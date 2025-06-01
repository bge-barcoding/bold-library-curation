#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use DBI;

# Simple loader for haplotype data into bold_haplotypes table
# Directly imports the TSV output from assess_haplotypes.pl

my $db_file;
my $tsv_file;
my $log_level = 'INFO';

GetOptions(
    'db=s'    => \$db_file,
    'tsv=s'   => \$tsv_file,
    'log=s'   => \$log_level,
);

die "Usage: $0 --db DATABASE --tsv HAPLOTYPE_TSV [--log LOG_LEVEL]\n" 
    unless $db_file && $tsv_file;

sub log_msg {
    my ($level, $msg) = @_;
    return if $log_level eq 'ERROR' && $level ne 'ERROR';
    return if $log_level eq 'WARN' && $level eq 'INFO';
    
    my $timestamp = scalar(localtime);
    print STDERR "[$timestamp] $level: $msg\n";
}

# Connect to database
log_msg('INFO', "Connecting to database: $db_file");
my $dbh = DBI->connect("dbi:SQLite:$db_file", "", "", {
    RaiseError => 1,
    AutoCommit => 1,  # We'll handle transactions manually
    sqlite_see_if_its_a_number => 1,
}) or die "Cannot connect to database: $DBI::errstr";

# Optimize SQLite for bulk loading
$dbh->do("PRAGMA journal_mode = WAL");
$dbh->do("PRAGMA synchronous = NORMAL"); 
$dbh->do("PRAGMA cache_size = -64000");  # 64MB cache
$dbh->do("PRAGMA temp_store = MEMORY");

# Clear existing haplotype data
log_msg('INFO', "Clearing existing haplotype data...");
$dbh->do("DELETE FROM bold_haplotypes");

# Prepare insert statement
my $insert_stmt = $dbh->prepare(q{
    INSERT OR REPLACE INTO bold_haplotypes (recordid, haplotype_id)
    VALUES (?, ?)
});

# Read and process haplotype TSV file
log_msg('INFO', "Processing haplotype file: $tsv_file");
open my $fh, '<', $tsv_file or die "Cannot open $tsv_file: $!";

my $header = <$fh>;  # Skip header line
chomp $header;
log_msg('INFO', "Header: $header");

my $processed_count = 0;
my $inserted_count = 0;
my $batch_size = 10000;  # Commit every 10k records for performance

# Start first transaction
$dbh->begin_work;

while (my $line = <$fh>) {
    chomp $line;
    next unless $line;
    
    my ($recordid, $haplotype_id, $status, $notes) = split /\t/, $line;
    
    # Skip failed assignments
    next unless $status && $status == 1;
    
    # Insert the haplotype assignment
    $insert_stmt->execute($recordid, $haplotype_id);
    
    $processed_count++;
    $inserted_count++;
    
    # Commit in batches for better performance with large datasets
    if ($inserted_count % $batch_size == 0) {
        $dbh->commit;
        log_msg('INFO', "Processed $processed_count records, inserted $inserted_count assignments...");
        $dbh->begin_work;  # Start next batch
    }
}

# Commit any remaining records
$dbh->commit;
close $fh;

# Generate statistics
my $total_assignments = $dbh->selectrow_array("SELECT COUNT(*) FROM bold_haplotypes");
my $unique_haplotypes = $dbh->selectrow_array("SELECT COUNT(DISTINCT haplotype_id) FROM bold_haplotypes");

log_msg('INFO', "Haplotype import completed successfully:");
log_msg('INFO', "  Total assignments imported: $total_assignments");
log_msg('INFO', "  Unique haplotypes: $unique_haplotypes");

# Show some example haplotype types
log_msg('INFO', "Sample haplotype assignments:");
my $samples = $dbh->selectall_arrayref(q{
    SELECT haplotype_id, COUNT(*) as count 
    FROM bold_haplotypes 
    GROUP BY haplotype_id 
    ORDER BY count DESC 
    LIMIT 10
}, { Slice => {} });

for my $sample (@$samples) {
    log_msg('INFO', "  $sample->{haplotype_id}: $sample->{count} records");
}

$dbh->disconnect;
log_msg('INFO', "Haplotype loading completed successfully");
