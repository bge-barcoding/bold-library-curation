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
if ($db_file) {
    $log->info("Connecting to database: $db_file");
    $io = BCDM::IO->new(db => $db_file);
} elsif ($tsv_file) {
    $log->info("Opening TSV file: $tsv_file");
    $io = BCDM::IO->new(tsv => $tsv_file);
}

# Process records in batches
$io->prepare_rs;
{
    my @queue;
    while (my $record = $io->next) {
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

