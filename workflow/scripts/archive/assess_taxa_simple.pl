#!/usr/bin/env perl
use strict;
use warnings;
use BCDM::IO;
use BCDM::ORM;
use Getopt::Long;
use Time::HiRes qw(time);

# Ultra-clean BAGS analysis with minimal, informative progress reporting
# Designed to provide clear visibility into processing without debug noise

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Command line options
my $db_file;
my $log_level = 'FATAL';  # Suppress BCDM debug output
my $progress_every = 50;  # Report every N species
my $quiet = 0;           # Minimal output mode

GetOptions(
    'db=s'        => \$db_file,
    'log=s'       => \$log_level,
    'progress=i'  => \$progress_every,
    'quiet'       => \$quiet,
);

die "Database file required (--db)\n" unless $db_file;

# Simple progress function
sub log_msg {
    my $msg = shift;
    return if $quiet;
    printf STDERR "[%s] %s\n", scalar(localtime), $msg;
}

# Force minimal logging from BCDM modules
BEGIN {
    $ENV{BCDM_LOG_LEVEL} = 'FATAL';
    # Redirect any remaining debug output to /dev/null equivalent
    if ($^O eq 'MSWin32') {
        $ENV{BCDM_LOG_FILE} = 'NUL';
    } else {
        $ENV{BCDM_LOG_FILE} = '/dev/null';
    }
}

use BCDM::BAGS;

# Connect with optimized settings
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

# Process species
my $iterator = $species_rs->search({}, { order_by => 'taxonid' });
while (my $taxon = $iterator->next) {
    # Create BAGS object (suppress any constructor debug output)
    my $bags;
    {
        local $SIG{__WARN__} = sub {}; # Suppress warnings during object creation
        $bags = BCDM::BAGS->new($taxon);
    }
    
    my $grade = $bags->grade;
    $grades{$grade}++;
    
    # Build result row
    my @lineage = reverse grep { $_->level =~ /^(order|family|genus)$/ } $taxon->lineage;
    my @row = (
        $taxon->taxonid,
        (map { $_->name } @lineage),
        $taxon->name,
        $grade
    );
    
    # Output BIN data
    for my $bin (@{ $bags->bins }) {
        next unless defined $bin && $bin =~ /^BOLD:/;
        my @sharers = $bags->taxa_sharing_bin($bin);
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
