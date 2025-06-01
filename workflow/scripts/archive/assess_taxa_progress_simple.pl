#!/usr/bin/env perl
use strict;
use warnings;
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;
use Getopt::Long;
use Time::HiRes qw(time);

# Simplified BAGS analysis with clean progress reporting
# Author: Improved version for clear progress tracking

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Process command line arguments
my $db_file;
my $log_level = 'FATAL';  # Suppress all debug/info output from BCDM modules
my $progress_interval = 25;  # Report progress every N species
my $milestone_interval = 100; # Major milestone every N species

GetOptions(
    'db=s'           => \$db_file,
    'log=s'          => \$log_level,
    'progress=i'     => \$progress_interval,
    'milestone=i'    => \$milestone_interval,
);

die "Database file required (--db)" unless $db_file;

# Clean progress logging to STDERR
sub progress_log {
    my $message = shift;
    my $timestamp = sprintf "[%s]", scalar localtime();
    print STDERR "$timestamp $message\n";
}

# Suppress the verbose BCDM logging by setting environment
$ENV{BCDM_LOG_LEVEL} = $log_level;

# Connect to database
progress_log("BAGS Analysis Starting");
progress_log("Database: $db_file");

my $orm = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { 
    quote_char => '"',
    on_connect_do => [
        'PRAGMA synchronous = NORMAL',
        'PRAGMA cache_size = -64000',
        'PRAGMA temp_store = MEMORY'
    ]
});

# Get total count for progress calculation
progress_log("Counting species to analyze...");
my $total_species = $orm->resultset('Taxa')->search({ 
    level => 'species', 
    name => { '!=' => '' } 
})->count;

progress_log("Found $total_species species to analyze");
progress_log("Progress will be reported every $progress_interval species");
progress_log("Major milestones every $milestone_interval species");
progress_log("Starting analysis...");
print STDERR "\n";

# Initialize tracking variables
my $processed = 0;
my $start_time = time();
my %grade_counts = ();
my $last_progress_time = $start_time;

# Print TSV header to STDOUT
my @header = qw[ taxonid order family genus species BAGS BIN sharers ];
print join("\t", @header), "\n";

# Get taxa iterator with proper ordering for predictable progress
my $taxa = $orm->resultset('Taxa')->search(
    { level => 'species', name => { '!=' => '' } },
    { order_by => 'name' }  # Alphabetical order for consistency
);

# Process each species
while (my $taxon = $taxa->next) {
    my $bags = BCDM::BAGS->new($taxon);
    my $grade = $bags->grade;
    $grade_counts{$grade}++;
    
    # Prepare basic result data
    my @result = (
        $taxon->taxonid, 
        map( { $_->name } ofg_lineage($taxon) ), 
        $taxon->name, 
        $grade
    );
    
    # Output BIN information
    BIN: for my $bin ( @{ $bags->bins } ) {
        next BIN unless defined $bin and $bin =~ /^BOLD:/;
        print join "\t", @result, $endpoint . $bin;
        my @sharers = $bags->taxa_sharing_bin($bin);
        print "\t", join(',', @sharers);
        print "\n";
    }
    
    $processed++;
    
    # Progress reporting
    if ($processed % $progress_interval == 0) {
        my $current_time = time();
        my $elapsed = $current_time - $start_time;
        my $rate = $processed / $elapsed;
        my $eta_seconds = ($total_species - $processed) / $rate;
        
        # Format ETA nicely
        my $eta_str = format_duration($eta_seconds);
        my $percent = sprintf("%.1f", ($processed / $total_species) * 100);
        my $rate_per_min = sprintf("%.1f", $rate * 60);
        
        progress_log(sprintf("Progress: %d/%d (%s%%) | Rate: %s/min | ETA: %s", 
                            $processed, $total_species, $percent, $rate_per_min, $eta_str));
        
        # Show grade distribution compactly
        my @grade_summary = map { "$_:$grade_counts{$_}" } sort keys %grade_counts;
        progress_log("Grades: " . join(", ", @grade_summary));
        
        $last_progress_time = $current_time;
    }
    
    # Major milestones
    if ($processed % $milestone_interval == 0) {
        my $elapsed_total = time() - $start_time;
        progress_log("*** MILESTONE: $processed species completed in " . 
                    format_duration($elapsed_total) . " ***");
        print STDERR "\n";  # Extra space for readability
    }
}

# Final summary
my $total_time = time() - $start_time;
my $final_rate = $processed / $total_time;

print STDERR "\n";
progress_log("=== BAGS ANALYSIS COMPLETED ===");
progress_log(sprintf("Total processed: %d species", $processed));
progress_log(sprintf("Total time: %s", format_duration($total_time)));
progress_log(sprintf("Average rate: %.1f species/minute", $final_rate * 60));
print STDERR "\n";

progress_log("FINAL GRADE DISTRIBUTION:");
my $total_graded = 0;
foreach my $grade (sort keys %grade_counts) {
    my $count = $grade_counts{$grade};
    my $percent = sprintf("%.1f", ($count / $processed) * 100);
    progress_log(sprintf("  Grade %s: %d species (%s%%)", $grade, $count, $percent));
    $total_graded += $count;
}
progress_log("Analysis complete!");

# Helper function to format duration nicely
sub format_duration {
    my $seconds = shift;
    return "0s" if $seconds < 1;
    
    my $hours = int($seconds / 3600);
    my $minutes = int(($seconds % 3600) / 60);
    my $secs = int($seconds % 60);
    
    if ($hours > 0) {
        return sprintf("%dh%02dm%02ds", $hours, $minutes, $secs);
    } elsif ($minutes > 0) {
        return sprintf("%dm%02ds", $minutes, $secs);
    } else {
        return sprintf("%ds", $secs);
    }
}

# Helper function for taxonomic lineage
sub ofg_lineage {
    my $taxon = shift;
    my @lineage = $taxon->lineage;
    my %keep = map { $_ => 1 } qw[ order family genus ];
    return grep { $keep{$_->level} } reverse @lineage;
}
