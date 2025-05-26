use strict;
use warnings;
use JSON;
use BCDM::IO;
use BCDM::ORM;
use BCDM::Criteria;
use LWP::UserAgent;
use Data::Dumper;
use Getopt::Long;
use Time::HiRes qw(usleep);
use Log::Log4perl qw(:easy);

# Adjustable sleep time (microseconds) and batch size
our $SLEEP = 500;
my $BATCH_SIZE = 1000;  # Reduced batch size to avoid 414 error

# API endpoints
my $base_url  = 'https://caos.boldsystems.org/api/images?processids=';
my $image_url = 'https://caos.boldsystems.org/api/objects/';

# Command-line arguments
my $db_file;
my $tsv_file;
my $log_level = 'INFO';
my $persist   = 0;
GetOptions(
    'db=s'         => \$db_file,
    'tsv=s'        => \$tsv_file,
    'log=s'        => \$log_level,
    'persist'      => \$persist,
    'batch_size=i' => \$BATCH_SIZE,  # Allow batch size to be configured
);

# Initialize Logger
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = $log_level, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# Instantiate logger
my $log = Log::Log4perl->get_logger('assess_criteria');
$log->info("Starting image retrieval with batch size: $BATCH_SIZE");

# Instantiate user agent
my $ua = LWP::UserAgent->new;
$log->info("UserAgent initialized");

# Open database or TSV
my $io;
my $schema; # For target filtering
if ($db_file) {
    $log->info("Connecting to database: $db_file");
    $io = BCDM::IO->new(db => $db_file);
    $schema = BCDM::ORM->connect("dbi:SQLite:$db_file");
} elsif ($tsv_file) {
    $log->info("Opening TSV file: $tsv_file");
    $io = BCDM::IO->new(tsv => $tsv_file);
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

# Helper function to check if record matches target taxa
sub record_matches_targets {
    my ($record) = @_;
    
    # If not using target filter, process all records
    return 1 unless $use_target_filter;
    
    # Check if record's taxon is in our target set
    my $taxonid = $record->taxonid;
    return 0 unless defined $taxonid;
    
    return exists $target_taxa{$taxonid};
}

# Process records in batches
my $total_processed = 0;
my $target_processed = 0;
$io->prepare_rs;
{
    my @queue;
    while (my $record = $io->next) {
        $total_processed++;
        
        # Skip records that don't match target taxa when target filtering is enabled
        unless ( record_matches_targets($record) ) {
            next;
        }
        
        $target_processed++;
        push @queue, $record;

        # Process batch when full
        if (scalar(@queue) >= $BATCH_SIZE) {
            eval { process_queue(@queue) };
            if ($@) {
                $log->error("Batch processing failed: $@");
                die $@;
            }

            # Clear queue and wait briefly
            @queue = ();
            usleep($SLEEP);
        }
    }

    # Process any remaining records
    process_queue(@queue) if @queue;
}

# Log summary statistics
if ( $use_target_filter ) {
    $log->info("Target filtering summary: processed $target_processed of $total_processed total records");
}
else {
    $log->info("Processed $total_processed records (no target filtering)");
}

sub process_queue {
    my @queue = @_;
    return unless @queue;  # Skip empty batches

    my @sub_batches;
    my $current_batch = '';

    # Split process IDs into smaller chunks that fit within a safe URL length
    foreach my $record (@queue) {
        my $pid = $record->processid;
        
        # Estimate if adding this PID would make the URL too long
        if (length($current_batch) + length($pid) + 1 > 7500) {  # 7500 chars as a safe limit
            push @sub_batches, $current_batch;
            $current_batch = $pid;
        } else {
            $current_batch .= ',' if $current_batch;
            $current_batch .= $pid;
        }
    }
    push @sub_batches, $current_batch if $current_batch;

    # Process each sub-batch separately
    foreach my $batch (@sub_batches) {
        my $wspoint = $base_url . $batch;
        $log->info("Fetching images for batch of up to " . scalar(@queue) . " records");

        # Send API request
        my $response = $ua->get($wspoint);
        
        # Handle API response
        if ($response->is_success) {
            my $json = $response->decoded_content;
            my $array_ref = decode_json $json;
            $log->debug(Dumper($array_ref));

            # Process each record in the batch
            for my $record (@queue) {
                my $pid = $record->processid;
                my $rid = $record->recordid;
                my ($match) = grep { $_->{processid} eq $pid } @$array_ref;
                my @result = ( $rid, $BCDM::Criteria::HAS_IMAGE );

                if ($match) {
                    push @result, 1, $image_url . $match->{objectid};
                } else {
                    push @result, 0, ':-(';
                }

                print join("\t", @result), "\n";
            }
        } else {
            $log->error("API request failed: " . $response->status_line);
            die $response->status_line;
        }
    }
}

