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
my $MAX_RETRIES = 5;     # Maximum retry attempts for API calls
my $RETRY_DELAY = 5;     # Initial retry delay in seconds (exponential backoff)
my $BATCH_DELAY = 1;     # Delay between batches in seconds (rate limiting)
my $CHECKPOINT_INTERVAL = 10000;  # Save checkpoint every N records

# API endpoints
my $base_url  = 'https://caos.boldsystems.org/api/images?processids=';
my $image_url = 'https://caos.boldsystems.org/api/objects/';

# Command-line arguments
my $db_file;
my $tsv_file;
my $output_file;
my $log_level = 'INFO';
my $persist   = 0;
my $resume    = 0;  # Flag to enable resuming from checkpoint
GetOptions(
    'db=s'         => \$db_file,
    'tsv=s'        => \$tsv_file,
    'output=s'     => \$output_file,
    'log=s'        => \$log_level,
    'persist'      => \$persist,
    'resume'       => \$resume,
    'batch_size=i' => \$BATCH_SIZE,     # Allow batch size to be configured
    'max_retries=i'=> \$MAX_RETRIES,    # Allow retry attempts to be configured
    'retry_delay=i'=> \$RETRY_DELAY,    # Allow retry delay to be configured
    'batch_delay=i'=> \$BATCH_DELAY,    # Allow batch delay to be configured
    'checkpoint_interval=i' => \$CHECKPOINT_INTERVAL,  # Allow checkpoint interval to be configured
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

# Setup output file handling
my $temp_output_file;
my $output_fh;

if ($output_file) {
    $temp_output_file = "$output_file.tmp";
    
    # Determine output mode based on resume flag and temp file existence
    my $output_mode = '>';  # Default: overwrite
    if ($resume && -f $temp_output_file) {
        $output_mode = '>>';  # Append mode when resuming
        $log->info("Resuming: appending to existing temporary output file: $temp_output_file");
    } else {
        $log->info("Writing to temporary output file: $temp_output_file");
    }
    
    open $output_fh, $output_mode, $temp_output_file 
        or die "Cannot open temporary output file $temp_output_file: $!";
    
    # Enable autoflush for immediate writes
    my $old_fh = select($output_fh);
    $| = 1;
    select($old_fh);
    
    # Write header only if creating new file
    if ($output_mode eq '>') {
        print $output_fh "recordid\tcriterion\tvalue\tnote\n";
    }
} else {
    # Fall back to stdout if no output file specified
    $output_fh = \*STDOUT;
    print $output_fh "recordid\tcriterion\tvalue\tnote\n";
    $log->info("Writing to stdout");
}

# Setup checkpoint file path
my $checkpoint_file;
my $checkpoint_data_file;
my $skip_count = 0;
my %processed_recordids;  # Track which records we've already processed

if ($db_file) {
    $checkpoint_file = "$db_file.image_checkpoint";
    $checkpoint_data_file = "$db_file.image_checkpoint.data";
    
    # Load checkpoint if resume flag is set
    if ($resume && -f $checkpoint_file) {
        open my $fh, '<', $checkpoint_file or die "Cannot open checkpoint file: $!";
        $skip_count = <$fh>;
        chomp $skip_count;
        close $fh;
        
        # Load processed record IDs
        if (-f $checkpoint_data_file) {
            open my $dfh, '<', $checkpoint_data_file or die "Cannot open checkpoint data: $!";
            while (my $line = <$dfh>) {
                chomp $line;
                $processed_recordids{$line} = 1;
            }
            close $dfh;
            $log->info("Loaded checkpoint: skipping first $skip_count records (" . 
                      (scalar keys %processed_recordids) . " already processed)");
        }
    } elsif ($resume && !-f $checkpoint_file) {
        $log->info("Resume flag set but no checkpoint found - starting from beginning");
    }
}

# Instantiate user agent with SSL configuration
my $ua = LWP::UserAgent->new(
    ssl_opts => {
        SSL_version => 'TLSv12:!SSLv2:!SSLv3',  # Force TLS 1.2+
        verify_hostname => 1,
        SSL_verify_mode => 0x00,  # More lenient verification
    },
    timeout => 300,  # 5 minute timeout
    agent => 'BOLD-Library-Curation/1.0',
);
$log->info("UserAgent initialized with enhanced SSL settings");

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
                # Use get_column to get the actual integer taxonid, not the Taxa object
                my $taxonid = $bt->get_column('taxonid');
                if (defined $taxonid) {
                    $target_taxa{$taxonid} = 1;
                    push @debug_target_list, $taxonid;
                    $bt_count++;
                    $log->debug("DEBUG: Loaded target taxonid: $taxonid");
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

# Helper function to check if record matches target taxa
sub record_matches_targets {
    my ($record) = @_;
    
    # If not using target filter, process all records
    return 1 unless $use_target_filter;
    
    # Get record info for debugging
    my $recordid = eval { $record->recordid } || "UNKNOWN";
    my $processid = eval { $record->processid } || "UNKNOWN";
    
    # Get the actual integer taxonid - use get_column if it's an ORM object
    my $taxonid;
    if (ref($record) && $record->can('get_column')) {
        $taxonid = eval { $record->get_column('taxonid') };
    } else {
        $taxonid = eval { $record->taxonid };
    }
    
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

# Helper function to save checkpoint
sub save_checkpoint {
    my ($count, @recordids) = @_;
    
    return unless $checkpoint_file;  # Only save if we have a checkpoint file
    
    # Save the count
    open my $fh, '>', $checkpoint_file or do {
        $log->warn("Cannot save checkpoint: $!");
        return;
    };
    print $fh $count;
    close $fh;
    
    # Append new record IDs to data file
    if (@recordids && $checkpoint_data_file) {
        open my $dfh, '>>', $checkpoint_data_file or do {
            $log->warn("Cannot save checkpoint data: $!");
            return;
        };
        print $dfh join("\n", @recordids), "\n";
        close $dfh;
    }
    
    $log->debug("Checkpoint saved: $count records processed");
}

# Process records in batches
my $total_processed = 0;
my $target_processed = 0;
my @processed_recordids = (); # Track processed record IDs for summary
my @checkpoint_batch = ();    # Batch of record IDs for checkpoint saving

$log->info("=== Starting image assessment record processing ===");
if ($skip_count > 0) {
    $log->info("Resuming from checkpoint - will skip first $skip_count records");
}

$io->prepare_rs;
{
    my @queue;
    while (my $record = $io->next) {
        $total_processed++;
        
        # Skip records if resuming from checkpoint
        if ($total_processed <= $skip_count) {
            next;
        }
        
        # Log progress periodically
        if ($total_processed % 1000 == 0) {
            $log->info("Progress: processed $total_processed total records, $target_processed targets matched");
        }
        
        # Check if record matches targets
        unless ( record_matches_targets($record) ) {
            next;
        }
        
        # Get record ID
        my $recordid = eval { $record->recordid } || "UNKNOWN";
        
        # Skip if already processed (from checkpoint)
        if (exists $processed_recordids{$recordid}) {
            $log->debug("Skipping already processed record: $recordid");
            next;
        }
        
        $target_processed++;
        my $processid = eval { $record->processid } || "UNKNOWN";
        my $taxonid;
        if (ref($record) && $record->can('get_column')) {
            $taxonid = eval { $record->get_column('taxonid') };
        } else {
            $taxonid = eval { $record->taxonid };
        }
        $taxonid = $taxonid || "UNDEFINED";
        
        push @processed_recordids, $recordid;
        push @checkpoint_batch, $recordid;
        $log->info("Processing target record for images #$target_processed: recordid=$recordid, processid=$processid, taxonid=$taxonid");
        
        push @queue, $record;

        # Process batch when full
        if (scalar(@queue) >= $BATCH_SIZE) {
            eval { process_queue(@queue) };
            if ($@) {
                $log->error("Batch processing failed: $@");
                # Save checkpoint before dying
                save_checkpoint($total_processed, @checkpoint_batch);
                die $@;
            }

            # Clear queue and wait briefly
            @queue = ();
            usleep($SLEEP);
            
            # Save checkpoint periodically
            if ($target_processed % $CHECKPOINT_INTERVAL == 0) {
                save_checkpoint($total_processed, @checkpoint_batch);
                @checkpoint_batch = ();  # Clear checkpoint batch after saving
                $log->info("Checkpoint saved at $total_processed records ($target_processed targets)");
            }
        }
    }

    # Process any remaining records
    if (@queue) {
        eval { process_queue(@queue) };
        if ($@) {
            $log->error("Final batch processing failed: $@");
            save_checkpoint($total_processed, @checkpoint_batch);
            die $@;
        }
    }
    
    # Save final checkpoint
    save_checkpoint($total_processed, @checkpoint_batch);
}

# Log summary statistics
$log->info("=== IMAGE ASSESSMENT SUMMARY ===");
if ( $use_target_filter ) {
    $log->info("Target filtering summary: processed $target_processed of $total_processed total records");
    $log->info("Target taxa used: [" . join(", ", @debug_target_list) . "]");
    $log->info("Processed record IDs: [" . join(", ", @processed_recordids) . "]");
    
    # Validate results
    if ($target_processed == 0) {
        $log->error("CRITICAL: No target records were processed for image assessment! Check target filtering logic.");
        $log->error("Expected target taxonids: [" . join(", ", @debug_target_list) . "]");
    } elsif ($target_processed < 10) {
        $log->warn("WARNING: Very few target records processed for images ($target_processed). Expected more based on database contents.");
    }
}
else {
    $log->info("Processed $total_processed records for image assessment (no target filtering)");
}

# Close output file handle
if ($output_file) {
    close $output_fh or $log->warn("Could not close output file handle: $!");
    
    # Rename temporary file to final output on success
    if (-f $temp_output_file) {
        rename $temp_output_file, $output_file 
            or die "Cannot rename temporary file to final output: $!";
        $log->info("Output file finalized: $output_file");
    }
}

# Clean up checkpoint files on successful completion
if ($checkpoint_file && -f $checkpoint_file) {
    unlink $checkpoint_file;
    $log->info("Checkpoint file removed after successful completion");
}
if ($checkpoint_data_file && -f $checkpoint_data_file) {
    unlink $checkpoint_data_file;
    $log->info("Checkpoint data file removed after successful completion");
}

$log->info("=== END IMAGE ASSESSMENT ===");

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

        # Retry logic with exponential backoff
        my $response;
        my $success = 0;
        
        for my $attempt (1..$MAX_RETRIES) {
            $response = $ua->get($wspoint);
            
            if ($response->is_success) {
                $success = 1;
                last;  # Success, exit retry loop
            }
            
            my $error = $response->status_line;
            $log->warn("API request failed (attempt $attempt/$MAX_RETRIES): $error");
            
            if ($attempt < $MAX_RETRIES) {
                my $delay = $RETRY_DELAY * (2 ** ($attempt - 1));  # Exponential backoff
                $log->info("Retrying in $delay seconds...");
                sleep($delay);
            } else {
                $log->error("Max retries exceeded. Final error: $error");
                die $error;
            }
        }
        
        # Handle API response
        if ($success) {
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

                print $output_fh join("\t", @result), "\n";
            }
        } else {
            $log->error("API request failed: " . $response->status_line);
            die $response->status_line;
        }
        
        # Rate limiting: delay between batches
        sleep($BATCH_DELAY) if $BATCH_DELAY > 0;
    }
}

