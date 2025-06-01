#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(time);
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;

# Simple BAGS analysis that filters output to show only progress
# Accepts that BCDM modules will produce debug output, but filters it

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
    
    my @lineage = reverse grep { $_->level =~ /^(order|family|genus)$/ } $taxon->lineage;
    my @row = (
        $taxon->taxonid,
        (map { $_->name } @lineage),
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
