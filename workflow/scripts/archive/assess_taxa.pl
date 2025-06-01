use strict;
use warnings;
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;
use Getopt::Long;
use Log::Log4perl qw(:easy);

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Process command line arguments
my $db_file;  # where to access database file
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
my $log = Log::Log4perl->get_logger('assess_taxa');
$log->info("Going to assess taxa for BAGS");

# Connect to database, prepare result set of species
$log->info("Going to connect to database $db_file");
my $orm = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { quote_char => '"' });

# NOTE: the query below now fetches ALL taxa in the normalized taxa table as produced by
# load_taxonomy.pl. This is because, for the BAGS grading, the basic currency is the taxon,
# not the BOLD record. Hence, in this script we iterate through the taxa, not the barcodes.
# Nevertheless, this particular query needs further refinement:
# - We only care to do this for species where the barcodes are top 3 in quality grading,
# - We only want to assess proper species, not provisional names or 'sp.' epithets,
# - For BGE we only care about those species that are in the target list.
# Meeting those requirements will involve a combination of things, such as getting the
# gap list cleaned up and linked to the database, pre-filtering the input data on curation
# quality, and refining the search predicates in the construct below:
my $taxa = $orm->resultset('Taxa')->search({ level => 'species', name => { '!=' => '' } });
$log->info("Will assess " . $taxa->count . " species");

# Print header
my @header = qw[ taxonid order family genus species BAGS BIN sharers ];
print join("\t", @header), "\n";

# Iterate over taxa
while (my $taxon = $taxa->next) {
    my $bags = BCDM::BAGS->new($taxon);
    my $grade = $bags->grade;
    my @result = ($taxon->taxonid, map( { $_->name } ofg_lineage($taxon) ), $taxon->name, $grade);
    BIN: for my $bin ( @{ $bags->bins } ) {
	next BIN unless defined $bin and $bin =~ /^BOLD:/;
        print join "\t", @result, $endpoint . $bin;
        my @sharers = $bags->taxa_sharing_bin($bin);
        print "\t", join(',', @sharers);
        print "\n";
    }
}

sub ofg_lineage {
    my $taxon = shift;
    my @lineage = $taxon->lineage;
    my %keep = map { $_ => 1 } qw[ order family genus ];
    return grep { $keep{$_->level} } reverse @lineage;
}

