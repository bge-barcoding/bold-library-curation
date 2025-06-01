#!/usr/bin/env perl
use strict;
use warnings;
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;
use Getopt::Long;
use Time::HiRes qw(time);

# Progress wrapper for BAGS analysis
# Provides meaningful progress updates without debug clutter

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Process command line arguments
my $db_file;
my $log_level = 'ERROR';  # Suppress debug output, only show errors and our progress
GetOptions(
    'db=s'  => \$db_file,
    'log=s' => \$log_level,
);

# Simple logging to STDERR for progress
sub log_progress {
    my $message = shift;
    print STDERR "[" . localtime() . "] $message\n";
}

# Connect to database
log_progress("Connecting to optimized database: $db_file");
my $orm = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { quote_char => '"' });

# Get species count for progress tracking
my $total_species = $orm->resultset('Taxa')->search({ level => 'species', name => { '!=' => '' } })->count;
log_progress("Found $total_species species to analyze");

# Initialize counters
my $processed = 0;
my $start_time = time();
my %grade_counts = ();

# Print header to STDOUT
my @header = qw[ taxonid order family genus species BAGS BIN sharers ];
print join("\t", @header), "\n";

# Get taxa iterator
my $taxa = $orm->resultset('Taxa')->search({ level => 'species', name => { '!=' => '' } });

# Process species with progress updates
while (my $taxon = $taxa->next) {
    my $bags = BCDM::BAGS->new($taxon);
    my $grade = $bags->grade;
    $grade_counts{$grade}++;
    
    my @result = ($taxon->taxonid, map( { $_->name } ofg_lineage($taxon) ), $taxon->name, $grade);
    
    BIN: for my $bin ( @{ $bags->bins } ) {
        next BIN unless defined $bin and $bin =~ /^BOLD:/;
        print join "\t", @result, $endpoint . $bin;
        my @sharers = $bags->taxa_sharing_bin($bin);
        print "\t", join(',', @sharers);
        print "\n";
    }
    
    $processed++;
    
    # Progress updates every 50 species
    if ($processed % 50 == 0) {
        my $elapsed = time() - $start_time;
        my $rate = $processed / $elapsed;
        my $eta_seconds = ($total_species - $processed) / $rate;
        my $eta_minutes = int($eta_seconds / 60);
        my $eta_hours = int($eta_minutes / 60);
        $eta_minutes = $eta_minutes % 60;
        
        my $percent = sprintf("%.1f", ($processed / $total_species) * 100);
        log_progress(sprintf("Progress: %d/%d (%s%%) - Rate: %.1f species/min - ETA: %02d:%02d", 
                            $processed, $total_species, $percent, $rate * 60, $eta_hours, $eta_minutes));
        
        # Show current grade distribution
        my $grade_summary = join(", ", map { "$_:$grade_counts{$_}" } sort keys %grade_counts);
        log_progress("Current grades: $grade_summary");
    }
    
    # Milestone updates
    if ($processed % 500 == 0) {
        log_progress("=== Milestone: $processed species completed ===");
    }
}

# Final summary
my $total_time = time() - $start_time;
my $final_rate = $processed / $total_time;

log_progress("=== BAGS Analysis Complete ===");
log_progress(sprintf("Processed: %d species in %.1f minutes", $processed, $total_time / 60));
log_progress(sprintf("Average rate: %.1f species/minute", $final_rate * 60));
log_progress("Final grade distribution:");
foreach my $grade (sort keys %grade_counts) {
    log_progress(sprintf("  Grade %s: %d species", $grade, $grade_counts{$grade}));
}

sub ofg_lineage {
    my $taxon = shift;
    my @lineage = $taxon->lineage;
    my %keep = map { $_ => 1 } qw[ order family genus ];
    return grep { $keep{$_->level} } reverse @lineage;
}
