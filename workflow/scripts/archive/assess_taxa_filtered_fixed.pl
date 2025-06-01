#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(time);
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;

# Fixed BAGS analysis that ensures proper column alignment
# by always outputting order, family, genus in correct positions

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

my $db_file;
my $progress_every = 50;
my $quiet = 0;

GetOptions(
    'db=s'        => \$db_file,
    'progress=i'  => \$progress_every,
    'quiet'       => \$quiet,
);

die "Database file required (--db)\n" unless $db_file;

# Progress to main STDERR (for log capture)
sub progress_msg {
    my $msg = shift;
    return if $quiet;
    printf STDERR "PROGRESS: [%s] %s\n", scalar(localtime), $msg;
}

# Connect to database
progress_msg("Connecting to database: $db_file");
my $orm = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { 
    quote_char => '"',
    on_connect_do => [
        'PRAGMA journal_mode = WAL',
        'PRAGMA synchronous = NORMAL', 
        'PRAGMA cache_size = -128000'
    ]
});

# Count species
progress_msg("Counting species...");
my $species_rs = $orm->resultset('Taxa')->search({ 
    level => 'species', 
    name => { '!=' => '' } 
});
my $total = $species_rs->count;
progress_msg("Found $total species to process");

my $count = 0;
my $start_time = time();
my %grades = ();

# Output header to STDOUT
print join("\t", qw[taxonid order family genus species BAGS BIN sharers]), "\n";

my $iterator = $species_rs->search({}, { order_by => 'taxonid' });
while (my $taxon = $iterator->next) {
    # BCDM modules will output debug info to STDERR, but we'll add our progress markers
    my $bags = BCDM::BAGS->new($taxon);
    my $grade = $bags->grade;
    $grades{$grade}++;
    
    # FIXED: Build taxonomy levels explicitly with fallbacks
    my @lineage = $taxon->lineage;
    my %lineage_by_level = ();
    
    # Extract each level into hash for reliable lookup
    for my $level_taxon (@lineage) {
        $lineage_by_level{$level_taxon->level} = $level_taxon->name;
    }
    
    # Ensure we always have order, family, genus in correct positions
    # Use empty string as fallback for missing levels
    my $order_name  = $lineage_by_level{'order'}  || '';
    my $family_name = $lineage_by_level{'family'} || '';
    my $genus_name  = $lineage_by_level{'genus'}  || '';
    
    my @row = (
        $taxon->taxonid,
        $order_name,
        $family_name,
        $genus_name,
        $taxon->name,
        $grade
    );
    
    for my $bin (@{ $bags->bins }) {
        next unless defined $bin && $bin =~ /^BOLD:/;
        my @sharers = $bags->taxa_sharing_bin($bin);
        print join("\t", @row, $endpoint . $bin, join(',', @sharers)), "\n";
    }
    
    $count++;
    
    # Clear progress updates with distinctive markers
    if ($count % $progress_every == 0) {
        my $elapsed = time() - $start_time;
        my $rate = $count / $elapsed;
        my $eta = ($total - $count) / $rate;
        my $pct = sprintf("%.1f", 100 * $count / $total);
        
        progress_msg(sprintf("SPECIES PROGRESS: %d/%d (%s%%) - %.1f/min - ETA: %s", 
                           $count, $total, $pct, $rate * 60, format_time($eta)));
        
        my $dist = join(", ", map { "$_:$grades{$_}" } sort keys %grades);
        progress_msg("GRADES: $dist") unless $quiet;
    }
    
    if ($count % 500 == 0) {
        progress_msg("MILESTONE: $count species completed");
    }
}

my $total_time = time() - $start_time;
progress_msg("ANALYSIS COMPLETE: $count species in " . format_time($total_time));
progress_msg("FINAL RATE: " . sprintf("%.1f", ($count / $total_time) * 60) . " species/minute");

unless ($quiet) {
    progress_msg("FINAL GRADES:");
    for my $grade (sort keys %grades) {
        my $pct = sprintf("%.1f", 100 * $grades{$grade} / $count);
        progress_msg(sprintf("  Grade %s: %d (%s%%)", $grade, $grades{$grade}, $pct));
    }
}

sub format_time {
    my $sec = shift;
    my $h = int($sec / 3600);
    my $m = int(($sec % 3600) / 60);
    my $s = int($sec % 60);
    
    return $h > 0 ? sprintf("%dh%02dm%02ds", $h, $m, $s) :
           $m > 0 ? sprintf("%dm%02ds", $m, $s) :
                    sprintf("%ds", $s);
}
