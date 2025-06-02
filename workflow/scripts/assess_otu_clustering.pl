#!/usr/bin/env perl

# assess_otu_clustering.pl
# Script to perform OTU clustering using VSEARCH
# Outputs recordid and OTU_ID pairs

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
    print STDERR "  --threads INT      Number of threads for VSEARCH (default: 8)\n";
    print STDERR "  --temp-dir DIR     Temporary directory (default: temp_otu)\n";
    print STDERR "  --help             Show this help\n";
    exit 1;
}

# VSEARCH expects threshold as decimal (0.0-1.0), not percentage
my $vsearch_threshold = $threshold;

# Logging function
sub log_message {
    my ($level, $message) = @_;
    my $timestamp = scalar(localtime());
    print STDERR "[$timestamp] $level: $message\n";
}

# Function to find vsearch binary
sub find_vsearch_binary {
    # List of possible locations for vsearch
    my @possible_paths = (
        'vsearch',                           # In PATH
        '/usr/bin/vsearch',                  # System install
        '/usr/local/bin/vsearch',            # Local install
        '$CONDA_PREFIX/bin/vsearch',         # Conda environment
        'which vsearch 2>/dev/null',        # Use which command
    );
    
    # Try to expand environment variables
    foreach my $path (@possible_paths) {
        if ($path =~ /\$(\w+)/) {
            my $env_var = $1;
            if (exists $ENV{$env_var}) {
                $path =~ s/\$$env_var/$ENV{$env_var}/g;
            }
        }
        
        # For the 'which' command, execute it
        if ($path =~ /^which/) {
            my $result = `$path`;
            chomp $result;
            if ($result && -x $result) {
                return $result;
            }
            next;
        }
        
        # Check if binary exists and is executable
        if (-x $path) {
            return $path;
        }
    }
    
    # If still not found, try conda-specific methods
    if (exists $ENV{CONDA_PREFIX}) {
        my $conda_vsearch = "$ENV{CONDA_PREFIX}/bin/vsearch";
        if (-x $conda_vsearch) {
            return $conda_vsearch;
        }
    }
    
    return undef;
}

log_message('INFO', "Starting OTU clustering analysis");
log_message('INFO', "Database: $db_file");
log_message('INFO', "Threshold: $threshold ($vsearch_threshold%)");
log_message('INFO', "Threads: $threads");
log_message('INFO', "Temp directory: $temp_dir");

# Create temporary directory
make_path($temp_dir);

# Connect to database
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", "", "", {
    RaiseError => 1,
    AutoCommit => 1
}) or die "Cannot connect to database: $DBI::errstr\n";

log_message('INFO', "Connected to database");

# Extract sequences with valid nucleotide data
log_message('INFO', "Extracting sequences for clustering");

my $fasta_file = "$temp_dir/sequences.fasta";
open(my $fasta_fh, '>', $fasta_file) or die "Cannot open $fasta_file: $!";

# Query to get sequences with nucleotide data
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
my $filtered_count = 0;
my @record_ids;

while (my ($recordid, $sequence) = $sth->fetchrow_array()) {
    # Clean sequence - remove whitespace and convert to uppercase
    $sequence =~ s/\s+//g;  # Remove whitespace
    $sequence = uc($sequence);
    
    # Trim leading and trailing gaps (alignment artifacts) but preserve internal gaps as N
    $sequence =~ s/^-+//;   # Remove leading gaps
    $sequence =~ s/-+$//;   # Remove trailing gaps
    $sequence =~ s/-/N/g;   # Replace any remaining internal gaps with N (unknown nucleotide)
    
    # Only skip sequences with invalid characters (not IUPAC nucleotide codes)
    # VSEARCH accepts all IUPAC codes: ACGTURYSWKMDBHVN
    if ($sequence =~ /[^ACGTURYSWKMDBHVN]/) {
        $filtered_count++;
        next;
    }
    
    # Skip very short sequences
    if (length($sequence) < 100) {
        $filtered_count++;
        next;
    }
    
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
if ($filtered_count > 0) {
    log_message('INFO', "Filtered out $filtered_count sequences (invalid characters or too short)");
}

if ($seq_count == 0) {
    log_message('ERROR', "No valid sequences found for clustering");
    # Output header only
    print "recordid\tOTU_ID\n";
    exit 0;
}

# Run VSEARCH clustering
log_message('INFO', "Running VSEARCH clustering");

my $clusters_file = "$temp_dir/clusters.uc";

# Try to find vsearch binary in multiple locations
my $vsearch_binary = find_vsearch_binary();
if (!$vsearch_binary) {
    die "VSEARCH binary not found. Please ensure VSEARCH is installed and in PATH.\n";
}

my $vsearch_cmd = "$vsearch_binary --cluster_fast $fasta_file --id $vsearch_threshold --uc $clusters_file --threads $threads --quiet";

log_message('INFO', "VSEARCH binary: $vsearch_binary");
log_message('INFO', "VSEARCH command: $vsearch_cmd");

my $vsearch_result = system($vsearch_cmd);
if ($vsearch_result != 0) {
    die "VSEARCH clustering failed with exit code: $vsearch_result\n";
}

log_message('INFO', "VSEARCH clustering completed");

# Parse clustering results
log_message('INFO', "Parsing clustering results");

my %otu_assignments;
my %otu_centroids;
my $otu_counter = 1;

open(my $uc_fh, '<', $clusters_file) or die "Cannot open $clusters_file: $!";

while (my $line = <$uc_fh>) {
    chomp $line;
    my @fields = split(/\t/, $line);
    
    # VSEARCH UC format:
    # Type H = hit (member of cluster)
    # Type C = centroid (cluster representative)
    # Type S = singleton (single member cluster)
    
    my $type = $fields[0];
    my $query_id = $fields[8];
    
    if ($type eq 'C' || $type eq 'S') {
        # This is a centroid/singleton - create new OTU
        my $otu_id = sprintf("OTU_%06d", $otu_counter);
        $otu_assignments{$query_id} = $otu_id;
        $otu_centroids{$query_id} = $otu_id;
        $otu_counter++;
    } elsif ($type eq 'H') {
        # This is a hit - assign to existing OTU
        my $target_id = $fields[9];
        if (exists $otu_centroids{$target_id}) {
            $otu_assignments{$query_id} = $otu_centroids{$target_id};
        }
    }
}

close($uc_fh);

my $total_otus = $otu_counter - 1;
log_message('INFO', "Created $total_otus OTUs from $seq_count sequences");

# Output results
log_message('INFO', "Generating output");

print "recordid\tOTU_ID\n";

# Output assignments for all processed sequences
foreach my $recordid (@record_ids) {
    my $otu_id = $otu_assignments{$recordid} || "UNASSIGNED";
    print "$recordid\t$otu_id\n";
}

# Clean up
$dbh->disconnect();
remove_tree($temp_dir);

log_message('INFO', "OTU clustering analysis completed");
log_message('INFO', "Output format: recordid, OTU_ID");

# Summary statistics
my %otu_counts;
foreach my $otu_id (values %otu_assignments) {
    $otu_counts{$otu_id}++;
}

log_message('INFO', "Summary statistics:");
log_message('INFO', "Total sequences processed: $seq_count");
log_message('INFO', "Total OTUs created: $total_otus");

my @singleton_otus = grep { $otu_counts{$_} == 1 } keys %otu_counts;
my $singleton_count = scalar(@singleton_otus);
log_message('INFO', "Singleton OTUs: $singleton_count");

my $largest_otu_size = 0;
my $largest_otu_id = '';
foreach my $otu_id (keys %otu_counts) {
    if ($otu_counts{$otu_id} > $largest_otu_size) {
        $largest_otu_size = $otu_counts{$otu_id};
        $largest_otu_id = $otu_id;
    }
}
log_message('INFO', "Largest OTU: $largest_otu_id with $largest_otu_size sequences");

exit 0;
