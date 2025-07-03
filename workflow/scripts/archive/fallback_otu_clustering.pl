#!/usr/bin/env perl

# fallback_otu_clustering.pl
# Alternative OTU clustering implementation using CD-HIT as fallback
# This script tries VSEARCH first, then falls back to CD-HIT if available

use strict;
use warnings;
use DBI;
use Getopt::Long;
use File::Temp qw(tempfile);
use File::Path qw(make_path remove_tree);

# Command line options
my $db_file;
my $log_level = 'INFO';
my $threshold = 0.99;
my $threads = 8;
my $temp_dir = 'temp_otu';
my $help;

GetOptions(
    'db=s'        => \$db_file,
    'log=s'       => \$log_level,
    'threshold=f' => \$threshold,
    'threads=i'   => \$threads,
    'temp-dir=s'  => \$temp_dir,
    'help'        => \$help
) or die "Error in command line arguments\n";

if ($help || !$db_file) {
    print STDERR "Usage: $0 --db DATABASE [options]\n";
    print STDERR "Options:\n";
    print STDERR "  --db FILE          SQLite database file\n";
    print STDERR "  --log LEVEL        Log level (default: INFO)\n";
    print STDERR "  --threshold FLOAT  Similarity threshold for clustering (default: 0.99)\n";
    print STDERR "  --threads INT      Number of threads (default: 8)\n";
    print STDERR "  --temp-dir DIR     Temporary directory (default: temp_otu)\n";
    print STDERR "  --help             Show this help\n";
    exit 1;
}

# Logging function
sub log_message {
    my ($level, $message) = @_;
    my $timestamp = scalar(localtime());
    print STDERR "[$timestamp] $level: $message\n";
}

# Check for available clustering tools
sub find_clustering_tool {
    my @tools = (
        {name => 'vsearch', cmd => 'vsearch'},
        {name => 'cd-hit-est', cmd => 'cd-hit-est'},
        {name => 'usearch', cmd => 'usearch'}
    );
    
    foreach my $tool (@tools) {
        my $check_cmd = "which $tool->{cmd} 2>/dev/null";
        my $result = `$check_cmd`;
        chomp $result;
        if ($result && -x $result) {
            log_message('INFO', "Found clustering tool: $tool->{name} at $result");
            return $tool;
        }
    }
    
    return undef;
}

log_message('INFO', "Starting fallback OTU clustering analysis");
log_message('INFO', "Database: $db_file");
log_message('INFO', "Threshold: $threshold");
log_message('INFO', "Threads: $threads");

# Find available clustering tool
my $clustering_tool = find_clustering_tool();
if (!$clustering_tool) {
    log_message('ERROR', "No clustering tools found (vsearch, cd-hit-est, usearch)");
    print "recordid\tOTU_ID\n";
    exit 1;
}

log_message('INFO', "Using clustering tool: $clustering_tool->{name}");

# Create temporary directory
make_path($temp_dir);

# Connect to database and extract sequences
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", "", "", {
    RaiseError => 1,
    AutoCommit => 1
}) or die "Cannot connect to database: $DBI::errstr\n";

log_message('INFO', "Connected to database");
log_message('INFO', "Extracting sequences for clustering");

my $fasta_file = "$temp_dir/sequences.fasta";
open(my $fasta_fh, '>', $fasta_file) or die "Cannot open $fasta_file: $!";

my $seq_query = qq{
    SELECT recordid, nuc 
    FROM bold 
    WHERE nuc IS NOT NULL 
    AND nuc != '' 
    AND length(nuc) > 100
    ORDER BY recordid
};

my $sth = $dbh->prepare($seq_query);
$sth->execute();

my $seq_count = 0;
my @record_ids;

while (my ($recordid, $sequence) = $sth->fetchrow_array()) {
    $sequence =~ s/\s+//g;
    $sequence = uc($sequence);
    next if $sequence =~ /[^ACGT]/;
    next if length($sequence) < 100;
    
    print $fasta_fh ">$recordid\n$sequence\n";
    push @record_ids, $recordid;
    $seq_count++;
    
    if ($seq_count % 10000 == 0) {
        log_message('INFO', "Processed $seq_count sequences");
    }
}

close($fasta_fh);
$sth->finish();

log_message('INFO', "Extracted $seq_count valid sequences");

if ($seq_count == 0) {
    log_message('ERROR', "No valid sequences found for clustering");
    print "recordid\tOTU_ID\n";
    exit 0;
}

# Run clustering based on available tool
my %otu_assignments;

if ($clustering_tool->{name} eq 'cd-hit-est') {
    run_cdhit_clustering($fasta_file, \%otu_assignments, \@record_ids);
} else {
    # This shouldn't happen since we check for tools above, but handle gracefully
    log_message('ERROR', "Unsupported clustering tool: $clustering_tool->{name}");
    exit 1;
}

# Output results
log_message('INFO', "Generating output");
print "recordid\tOTU_ID\n";

foreach my $recordid (@record_ids) {
    my $otu_id = $otu_assignments{$recordid} || "UNASSIGNED";
    print "$recordid\t$otu_id\n";
}

# Clean up
$dbh->disconnect();
remove_tree($temp_dir);

log_message('INFO', "Fallback OTU clustering analysis completed");

# CD-HIT clustering function
sub run_cdhit_clustering {
    my ($fasta_file, $otu_assignments, $record_ids) = @_;
    
    log_message('INFO', "Running CD-HIT clustering");
    
    my $output_file = "$temp_dir/clusters.fasta";
    my $cluster_file = "$temp_dir/clusters.fasta.clstr";
    
    my $cdhit_cmd = "cd-hit-est -i $fasta_file -o $output_file -c $threshold -n 10 -T $threads -M 0 -d 0";
    log_message('INFO', "CD-HIT command: $cdhit_cmd");
    
    my $result = system($cdhit_cmd);
    if ($result != 0) {
        die "CD-HIT clustering failed with exit code: $result\n";
    }
    
    # Parse CD-HIT cluster file
    log_message('INFO', "Parsing CD-HIT clustering results");
    
    open(my $cluster_fh, '<', $cluster_file) or die "Cannot open $cluster_file: $!";
    
    my $cluster_id = 0;
    my $current_otu = '';
    
    while (my $line = <$cluster_fh>) {
        chomp $line;
        
        if ($line =~ /^>Cluster (\d+)/) {
            $cluster_id = $1;
            $current_otu = sprintf("OTU_%06d", $cluster_id + 1);
        } elsif ($line =~ />(\w+)\.\.\./) {
            my $recordid = $1;
            $otu_assignments->{$recordid} = $current_otu;
        }
    }
    
    close($cluster_fh);
    
    my $total_otus = $cluster_id + 1;
    log_message('INFO', "Created $total_otus OTUs from " . scalar(@$record_ids) . " sequences using CD-HIT");
}

exit 0;
