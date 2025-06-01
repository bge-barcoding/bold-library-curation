use strict;
use warnings;
use BCDM::ORM;
use Getopt::Long;
use Log::Log4perl qw(:easy);

# Levels in the BOLD taxonomy
my @levels = qw[
    kingdom
    phylum
    class
    order
    family
    subfamily
    genus
    species
    subspecies
];

# Process command line arguments
my $db_file; # where to access database file
my $log_level = 'INFO'; # verbosity level for logger
GetOptions(
    'db=s'  => \$db_file,
    'log=s' => \$log_level,
);

# Initialize Log::Log4perl
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = $log_level, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# Instantiate logger
my $log = Log::Log4perl->get_logger('load_taxonomy');

# Connect to the database. Because there is a column `order` in the BCDM, we need to
# escape this. In SQLite that's with double quotes.
$log->info("Going to connect to database $db_file");
my $schema = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { quote_char => '"' });

# Iterate over all BOLD records to link to the taxa table
$log->info("Going to link BOLD records to normalized taxa");
my $rs = $schema->resultset('Bold')->search({});
my $tree = {}; # this is a cache that is gradually filled up
while (my $record = $rs->next) {
    $log->info('Processing record '.$record->recordid) unless $record->recordid % 1000;

    # Find lowest defined taxon name in BOLD record by going backward
    my $max = get_lowest_defined_taxon($record);

    # Find lowest defined taxon object either from cache or DB
    my $taxonid = get_taxonid($record, $max, $tree);

    # Update the taxon id
    $record->update({ 'taxonid' => $taxonid });
}

sub get_taxonid {
    my ( $record, $max, $node ) = @_;
    my $taxonid = undef;
    my $kingdom = $record->kingdom;
    for my $i ( 0 .. $max ) {
        my $level = $levels[$i];
        my $name = $record->$level;

        # already seen this taxonomic path in the tree
        if ( $node->{$name} ) {
            $taxonid = $node->{$name}->{object};
        }

        # not seen yet
        else {

            # create node
            $taxonid = get_node($kingdom, $name, $level, $taxonid);

            # extend tree
            $node->{$name} = { children => {}, object => $taxonid };
        }

        # traverse down the tree
        $node = $node->{$name}->{children};
    }
    return $taxonid;
}

sub get_lowest_defined_taxon {
    my $record = shift;
    for ( my $i = $#levels; $i >= 0; $i-- ) {
        my $level = $levels[$i];
        if ( $record->$level ne 'None' ) {
            return $i;
        }
    }
}

sub get_node {
    my ( $kingdom, $name, $level, $parent ) = @_;

    # Compose search predicates
    my $predicates = { 'kingdom' => $kingdom, 'name' => $name, 'level' => $level };
    $predicates->{parent_taxonid} = $parent if defined $parent;

    # Either an object or undef
    if ( my $taxon = $schema->resultset('Taxa')->find($predicates) ) {
        return $taxon;
    }
    else {
        $taxon = $schema->resultset('Taxa')->create($predicates);
        $taxon->insert;
        return $taxon->taxonid;
    }
}

