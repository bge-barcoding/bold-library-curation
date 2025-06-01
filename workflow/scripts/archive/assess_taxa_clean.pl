#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(time);

# Ultra-clean BAGS analysis - completely suppresses BCDM debug output
# Uses output redirection to ensure no debug messages leak through

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Command line options
my $db_file;
my $progress_every = 50;
my $quiet = 0;

GetOptions(
    'db=s'        => \$db_file,
    'progress=i'  => \$progress_every,
    'quiet'       => \$quiet,
);

die "Database file required (--db)\n" unless $db_file;

sub log_msg {
    my $msg = shift;
    return if $quiet;
    printf STDERR "[%s] %s\n", scalar(localtime), $msg;
}

# Completely override Log4perl before any modules load
BEGIN {
    # Create a minimal Log4perl configuration that sends everything to /dev/null
    my $null_device = ($^O eq 'MSWin32') ? 'NUL' : '/dev/null';
    
    $ENV{LOG4PERL_CONF} = qq{
        log4perl.rootLogger = OFF
        log4perl.appender.Null = Log::Log4perl::Appender::File
        log4perl.appender.Null.filename = $null_device
        log4perl.appender.Null.layout = Log::Log4perl::Layout::SimpleLayout
    };
}

# Load Log4perl first with our configuration
use Log::Log4perl qw(:easy);
my $null_device = ($^O eq 'MSWin32') ? 'NUL' : '/dev/null';
Log::Log4perl->init(\qq{
    log4perl.rootLogger = OFF
    log4perl.appender.Null = Log::Log4perl::Appender::File
    log4perl.appender.Null.filename = $null_device
    log4perl.appender.Null.layout = Log::Log4perl::Layout::SimpleLayout
});

# Now load BCDM modules
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;

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

log_msg("Counting species...");
my $species_rs = $orm->resultset('Taxa')->search({ 
    level => 'species', 
    name => { '!=' => '' } 
});
my $total = $species_rs->count;
log_msg("Found $total species to process");

my $count = 0;
my $start_time = time();
my %grades = ();

# Output header
print join("\t", qw[taxonid order family genus species BAGS BIN sharers]), "\n";

my $iterator = $species_rs->search({}, { order_by => 'taxonid' });
while (my $taxon = $iterator->next) {
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
    
    if ($count % $progress_every == 0) {
        my $elapsed = time() - $start_time;
        my $rate = $count / $elapsed;
        my $eta = ($total - $count) / $rate;
        my $pct = sprintf("%.1f", 100 * $count / $total);
        
        log_msg(sprintf("Progress: %d/%d (%s%%) - %.1f/min - ETA: %s", 
                       $count, $total, $pct, $rate * 60, format_time($eta)));
        
        my $dist = join(", ", map { "$_:$grades{$_}" } sort keys %grades);
        log_msg("Grades: $dist") unless $quiet;
    }
    
    if ($count % 500 == 0) {
        log_msg("=== MILESTONE: $count species completed ===");
    }
}

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

sub format_time {
    my $sec = shift;
    my $h = int($sec / 3600);
    my $m = int(($sec % 3600) / 60);
    my $s = int($sec % 60);
    
    return $h > 0 ? sprintf("%dh%02dm%02ds", $h, $m, $s) :
           $m > 0 ? sprintf("%dm%02ds", $m, $s) :
                    sprintf("%ds", $s);
}
