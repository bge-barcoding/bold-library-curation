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
my @debug_target_list; # For debugging

if ( $db_file && $schema ) {
    eval {
        # Check if targets table exists and has data
        my $target_count = $schema->resultset('Target')->count;
        $log->info("DEBUG: Found $target_count target entries in database");
        
        if ( $target_count > 0 ) {
            $log->info("Found $target_count targets - enabling target filtering");
            $use_target_filter = 1;
            
            # Build set of target taxon IDs with explicit ordering and verification
            my $bold_targets = $schema->resultset('BoldTarget')->search({}, { order_by => 'taxonid' });
            my $bt_count = 0;
            
            while ( my $bt = $bold_targets->next ) {
                my $taxonid = $bt->taxonid;
                if (defined $taxonid) {
                    $target_taxa{$taxonid} = 1;
                    push @debug_target_list, $taxonid;
                    $bt_count++;
                } else {
                    $log->warn("WARNING: Found BoldTarget record with undefined taxonid");
                }
            }
            
            # Sort for consistent logging
            @debug_target_list = sort { $a <=> $b } @debug_target_list;
            
            my $hash_size = keys %target_taxa;
            $log->info("Built target filter with $hash_size target taxa");
            $log->info("DEBUG: Target taxonids: [" . join(", ", @debug_target_list) . "]");
            $log->info("DEBUG: BoldTarget records processed: $bt_count");
            
            # Verify consistency
            if ($hash_size != $bt_count) {
                $log->error("CRITICAL: Hash size ($hash_size) != BoldTarget count ($bt_count)");
            }
            if ($hash_size == 0) {
                $log->error("CRITICAL: No target taxa loaded - disabling target filtering");
                $use_target_filter = 0;
            }
        }
        else {
            $log->info("No targets found - processing all records");
        }
    };
    if ( $@ ) {
        $log->warn("Could not check for targets (tables may not exist): $@");
        $log->info("Processing all records");
        $use_target_filter = 0;
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
    
    # Get record info for debugging
    my $recordid = eval { $record->recordid } || "UNKNOWN";
    my $processid = eval { $record->processid } || "UNKNOWN";
    
    # Check if record's taxon is in our target set
    my $taxonid = eval { $record->taxonid };
    
    # Handle undefined taxonid
    if (!defined $taxonid) {
        $log->debug("DEBUG: Record $recordid (processid: $processid) has undefined taxonid - skipping");
        return 0;
    }
    
    # Check if this taxonid is in our target set
    my $is_target = exists $target_taxa{$taxonid};
    
    if ($is_target) {
        $log->debug("DEBUG: ✓ PROCESSING target record $recordid (processid: $processid, taxonid: $taxonid)");
        return 1;
    } else {
        $log->debug("DEBUG: ✗ SKIPPING non-target record $recordid (processid: $processid, taxonid: $taxonid)");
        return 0;
    }
}

# Iterate over all BOLD records
my $total_processed = 0;
my $target_processed = 0;
my @processed_recordids = (); # Track processed record IDs for summary

$log->info("=== Starting record processing ===");
$io->prepare_rs;
while (my $record = $io->next) {
    $total_processed++;
    
    # Log progress periodically
    if ($total_processed % 1000 == 0) {
        $log->info("Progress: processed $total_processed total records, $target_processed targets matched");
    }
    
    # Check if record matches targets
    unless ( record_matches_targets($record) ) {
        next;
    }
    
    $target_processed++;
    my $recordid = eval { $record->recordid } || "UNKNOWN";
    my $processid = eval { $record->processid } || "UNKNOWN";
    my $taxonid = eval { $record->taxonid } || "UNDEFINED";
    
    push @processed_recordids, $recordid;
    
    $log->info("Processing target record #$target_processed: recordid=$recordid, processid=$processid, taxonid=$taxonid");

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
$log->info("=== ASSESSMENT SUMMARY ===");
if ( $use_target_filter ) {
    $log->info("Target filtering summary: processed $target_processed of $total_processed total records");
    $log->info("Target taxa used: [" . join(", ", @debug_target_list) . "]");
    $log->info("Processed record IDs: [" . join(", ", @processed_recordids) . "]");
    
    # Validate results
    if ($target_processed == 0) {
        $log->error("CRITICAL: No target records were processed! Check target filtering logic.");
        $log->error("Expected target taxonids: [" . join(", ", @debug_target_list) . "]");
    } elsif ($target_processed < 10) {
        $log->warn("WARNING: Very few target records processed ($target_processed). Expected more based on database contents.");
    }
}
else {
    $log->info("Processed $total_processed records (no target filtering)");
}
$log->info("=== END ASSESSMENT ===");
