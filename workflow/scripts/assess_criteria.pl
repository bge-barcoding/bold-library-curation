use strict;
use warnings;
use BCDM::IO;
use BCDM::ORM;
use BCDM::Criteria;
use Getopt::Long;
use Log::Log4perl qw(:easy);

# Process command line arguments
my $db_file;  # where to access database file
my $tsv_file; # using the raw TSV instead
my $log_level = 'INFO'; # verbosity level for logger
my $persist   = 0;
my @criteria; # e.g. HAS_IMAGE
GetOptions(
    'db=s'       => \$db_file,
    'tsv=s'      => \$tsv_file,
    'log=s'      => \$log_level,
    'criteria=s' => \@criteria,
    'persist'    => \$persist,
);

# Initialize Log::Log4perl
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = $log_level, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# Instantiate logger
my $log = Log::Log4perl->get_logger('assess_criteria');
$log->info("Going to assess criteria: @criteria");

# Connect to the database. Because there is a column `order` in the BCDM, we need to
# escape this. In SQLite that's with double quotes.
my $io;
my $schema; # For target filtering
if ( $db_file ) {
    $log->info("Going to connect to database $db_file");
    $io = BCDM::IO->new( db => $db_file );
    $schema = BCDM::ORM->connect("dbi:SQLite:$db_file");
}
elsif ( $tsv_file) {
    $log->info("Going to open TSV file $tsv_file");
    $io = BCDM::IO->new( tsv => $tsv_file );
}

# Check if target list is being used and build target taxa set
my %target_taxa;
my $use_target_filter = 0;
if ( $db_file && $schema ) {
    eval {
        # Check if targets table exists and has data
        my $target_count = $schema->resultset('Target')->count;
        if ( $target_count > 0 ) {
            $log->info("Found $target_count targets - enabling target filtering");
            $use_target_filter = 1;
            
            # Build set of target taxon IDs
            my $bold_targets = $schema->resultset('BoldTarget');
            while ( my $bt = $bold_targets->next ) {
                $target_taxa{$bt->taxonid} = 1;
            }
            $log->info("Built target filter with " . keys(%target_taxa) . " target taxa");
        }
        else {
            $log->info("No targets found - processing all records");
        }
    };
    if ( $@ ) {
        $log->warn("Could not check for targets (tables may not exist): $@");
        $log->info("Processing all records");
    }
}

# Create map of loaded criteria
my %crit;
for my $c ( @criteria ) {
    next unless $c;
    $log->info("Attempting to load '$c'");
    my $impl = BCDM::Criteria->load_criterion($c);
    $crit{$c} = $impl;
    $log->info("Loaded implementation $c => $impl");
}

# Helper function to check if record matches target taxa
sub record_matches_targets {
    my ($record) = @_;
    
    # If not using target filter, process all records
    return 1 unless $use_target_filter;
    
    # Check if record's taxon is in our target set
    # The record should have taxonid field based on the database schema
    my $taxonid = $record->taxonid;
    return 0 unless defined $taxonid;
    
    return exists $target_taxa{$taxonid};
}

# Iterate over all BOLD records
my $total_processed = 0;
my $target_processed = 0;
$io->prepare_rs;
while (my $record = $io->next) {
    $total_processed++;
    $log->info("Processing record ".$record->recordid) unless $record->recordid % 10_000;
    
    # Skip records that don't match target taxa when target filtering is enabled
    unless ( record_matches_targets($record) ) {
        next;
    }
    
    $target_processed++;

    # Iterate over loaded modules
    for my $impl ( values %crit ) {

        # Code reference passed into the assess() function to run batches
        my $handler = sub {
            my ( $status, $notes, $i ) = @_;

            # Persist to database or print to stdout
            if ( $persist ) {
                $impl->persist(
                    record => $record,
                    status => $status,
                    notes  => $notes
                );
            }
            else {

                # No primary key in the output, needs to be generated
                my $cid = $impl->_criterion;
                my $rid = ( $record->recordid - $impl->_batch_size + $i );
                print join( "\t", $rid, $cid, $status, $notes ), "\n";
            }
        };

        # Do the assessment
        $impl->assess(
            record  => $record,
            handler => $handler
        );


    }
}

# Log summary statistics
if ( $use_target_filter ) {
    $log->info("Target filtering summary: processed $target_processed of $total_processed total records");
}
else {
    $log->info("Processed $total_processed records (no target filtering)");
}
