#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(time);

# Clean BAGS analysis with logging suppression
# This version properly suppresses BCDM logging before module loading

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Command line options
my $db_file;
my $log_level = 'FATAL';
my $progress_every = 50;
my $quiet = 0;

GetOptions(
    'db=s'        => \$db_file,
    'log=s'       => \$log_level,
    'progress=i'  => \$progress_every,
    'quiet'       => \$quiet,
);

die "Database file required (--db)\n" unless $db_file;

# Progress logging function
sub log_msg {
    my $msg = shift;
    return if $quiet;
    printf STDERR "[%s] %s\n", scalar(localtime), $msg;
}

# CRITICAL: Override Log::Log4perl configuration BEFORE loading BCDM modules
BEGIN {
    # Set up minimal logging configuration to suppress debug output
    $ENV{LOG4PERL_CONF} = q{
        log4perl.rootLogger = FATAL, Screen
        log4perl.appender.Screen = Log::Log4perl::Appender::Screen
        log4perl.appender.Screen.stderr = 1
        log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
    };
    
    # Also try alternative environment variables
    $ENV{BCDM_LOG_LEVEL} = 'FATAL';
    $ENV{LOG_LEVEL} = 'FATAL';
}

# Now load BCDM modules with suppressed logging
use Log::Log4perl qw(:easy);
Log::Log4perl->init(\<<"END");
    log4perl.rootLogger = FATAL, Screen
    log4perl.appender.Screen = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr = 1
    log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# Redirect STDERR temporarily during module loading to suppress any remaining output
{
    open my $old_stderr, '>&', \*STDERR or die "Can't dup STDERR: $!";
    close STDERR;
    open STDERR, '>', ($^O eq 'MSWin32' ? 'NUL' : '/dev/null') or die "Can't redirect STDERR: $!";
    
    # Load BCDM modules
    require BCDM::IO;
    require BCDM::ORM;
    require BCDM::BAGS;
    
    # Restore STDERR
    close STDERR;
    open STDERR, '>&', $old_stderr or die "Can't restore STDERR: $!";
}

# Connect to database
log_msg("Connecting to database: $db_file");
my $orm = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { 
    quote_char => '"',
    on_connect_do => [
        'PRAGMA journal_mode = WAL',
        'PRAGMA synchronous = NORMAL', 
        'PRAGMA cache_size = -128000',
        'PRAGMA temp_store = MEMORY'
    ]
});

# Count species
log_msg("Counting species...");
my $species_rs = $orm->resultset('Taxa')->search({ 
    level => 'species', 
    name => { '!=' => '' } 
});
my $total = $species_rs->count;
log_msg("Found $total species to process");

# Initialize tracking
my $count = 0;
my $start_time = time();
my %grades = ();

# Output header
print join("\t", qw[taxonid order family genus species BAGS BIN sharers]), "\n";

# Process species with output suppression
my $iterator = $species_rs->search({}, { order_by => 'taxonid' });
while (my $taxon = $iterator->next) {
    # Suppress all output during BAGS processing
    my ($bags, $grade, @bins);
    {
        # Temporarily redirect STDERR to suppress BCDM debug output
        open my $old_stderr, '>&', \*STDERR or die "Can't dup STDERR: $!";
        close STDERR;
        open STDERR, '>', ($^O eq 'MSWin32' ? 'NUL' : '/dev/null') or die "Can't redirect STDERR: $!";
        
        # Create BAGS object and get results
        $bags = BCDM::BAGS->new($taxon);
        $grade = $bags->grade;
        @bins = @{ $bags->bins };
        
        # Restore STDERR
        close STDERR;
        open STDERR, '>&', $old_stderr or die "Can't restore STDERR: $!";
    }
    
    $grades{$grade}++;
    
    # Build result row
    my @lineage = reverse grep { $_->level =~ /^(order|family|genus)$/ } $taxon->lineage;
    my @row = (
        $taxon->taxonid,
        (map { $_->name } @lineage),
        $taxon->name,
        $grade
    );
    
    # Output BIN data (suppress sharing calculation output too)
    for my $bin (@bins) {
        next unless defined $bin && $bin =~ /^BOLD:/;
        
        my @sharers;
        {
            # Suppress output during sharing calculation
            open my $old_stderr, '>&', \*STDERR or die "Can't dup STDERR: $!";
            close STDERR;
            open STDERR, '>', ($^O eq 'MSWin32' ? 'NUL' : '/dev/null') or die "Can't redirect STDERR: $!";
            
            @sharers = $bags->taxa_sharing_bin($bin);
            
            # Restore STDERR
            close STDERR;
            open STDERR, '>&', $old_stderr or die "Can't restore STDERR: $!";
        }
        
        print join("\t", @row, $endpoint . $bin, join(',', @sharers)), "\n";
    }
    
    $count++;
    
    # Progress updates
    if ($count % $progress_every == 0) {
        my $elapsed = time() - $start_time;
        my $rate = $count / $elapsed;
        my $eta = ($total - $count) / $rate;
        my $pct = sprintf("%.1f", 100 * $count / $total);
        
        log_msg(sprintf("Progress: %d/%d (%s%%) - %.1f/min - ETA: %s", 
                       $count, $total, $pct, $rate * 60, format_time($eta)));
        
        # Show grade distribution
        my $dist = join(", ", map { "$_:$grades{$_}" } sort keys %grades);
        log_msg("Grades: $dist") unless $quiet;
    }
    
    # Milestone markers
    if ($count % 500 == 0) {
        log_msg("=== MILESTONE: $count species completed ===");
    }
}

# Final summary
my $total_time = time() - $start_time;
log_msg("COMPLETED: $count species in " . format_time($total_time));
log_msg("Average rate: " . sprintf("%.1f", ($count / $total_time) * 60) . " species/minute");

unless ($quiet) {
    log_msg("Grade distribution:");
    for my $grade (sort keys %grades) {
        my $pct = sprintf("%.1f", 100 * $grades{$grade} / $count);
        log_msg(sprintf("  %s: %d (%s%%)", $grade, $grades{$grade}, $pct));
    }
}

# Helper function to format time duration
sub format_time {
    my $sec = shift;
    my $h = int($sec / 3600);
    my $m = int(($sec % 3600) / 60);
    my $s = int($sec % 60);
    
    return $h > 0 ? sprintf("%dh%02dm%02ds", $h, $m, $s) :
           $m > 0 ? sprintf("%dm%02ds", $m, $s) :
                    sprintf("%ds", $s);
}
