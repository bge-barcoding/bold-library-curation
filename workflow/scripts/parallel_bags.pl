#!/usr/bin/env perl
use strict;
use warnings;
use BCDM::IO;
use BCDM::ORM;
use BCDM::BAGS;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use POSIX qw(ceil);

my $endpoint = 'https://boldsystems.org/index.php/Public_BarcodeCluster?clusteruri=';

# Process command line arguments
my $db_file;
my $log_level = 'INFO';
my $chunk_id = 0;      # Which chunk to process (from SLURM_ARRAY_TASK_ID)
my $total_chunks = 1;  # Total number of chunks
my $output_dir = 'results/chunks';

GetOptions(
    'db=s'          => \$db_file,
    'log=s'         => \$log_level,
    'chunk-id=i'    => \$chunk_id,
    'total-chunks=i'=> \$total_chunks,
    'output-dir=s'  => \$output_dir,
);

# Use SLURM environment variables if available
$chunk_id = $ENV{SLURM_ARRAY_TASK_ID} if defined $ENV{SLURM_ARRAY_TASK_ID};
$total_chunks = $ENV{SLURM_ARRAY_TASK_MAX} - $ENV{SLURM_ARRAY_TASK_MIN} + 1 
    if defined $ENV{SLURM_ARRAY_TASK_MAX} && defined $ENV{SLURM_ARRAY_TASK_MIN};

# Initialize logging
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = $log_level, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

my $log = Log::Log4perl->get_logger("assess_taxa_chunk_$chunk_id");
$log->info("Processing chunk $chunk_id of $total_chunks");

# Connect to database
my $orm = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { 
    quote_char => '"',
    sqlite_busy_timeout => 30000,  # 30 second timeout for database locks
});

# Get total species count and calculate chunk boundaries
my $total_species = $orm->resultset('Taxa')->search({ 
    level => 'species', 
    name => { '!=' => '' } 
})->count;

my $chunk_size = ceil($total_species / $total_chunks);
my $offset = $chunk_id * $chunk_size;

$log->info("Total species: $total_species, chunk size: $chunk_size, offset: $offset");

# Get this chunk's species with LIMIT/OFFSET
my $taxa = $orm->resultset('Taxa')->search(
    { level => 'species', name => { '!=' => '' } },
    { 
        rows => $chunk_size,
        offset => $offset,
        order_by => 'taxonid'  # Ensure consistent ordering
    }
);

my $species_count = $taxa->count;
$log->info("Processing $species_count species in chunk $chunk_id");

# Create output directory
system("mkdir -p $output_dir") unless -d $output_dir;

# Open output file for this chunk
my $output_file = "$output_dir/bags_chunk_$chunk_id.tsv";
open my $fh, '>', $output_file or die "Cannot open $output_file: $!";

# Print header only for first chunk
if ($chunk_id == 0) {
    my @header = qw[ taxonid order family genus species BAGS BIN sharers ];
    print $fh join("\t", @header), "\n";
}

# Process species in this chunk
my $processed = 0;
while (my $taxon = $taxa->next) {
    $processed++;
    $log->info("Processing species $processed/$species_count: " . $taxon->name) if $processed % 100 == 0;
    
    my $bags = BCDM::BAGS->new($taxon);
    my $grade = $bags->grade;
    my @result = ($taxon->taxonid, map( { $_->name } ofg_lineage($taxon) ), $taxon->name, $grade);
    
    BIN: for my $bin ( @{ $bags->bins } ) {
        next BIN unless defined $bin and $bin =~ /^BOLD:/;
        print $fh join "\t", @result, $endpoint . $bin;
        my @sharers = $bags->taxa_sharing_bin($bin);
        print $fh "\t", join(',', @sharers);
        print $fh "\n";
    }
}

close $fh;
$log->info("Chunk $chunk_id completed. Output written to $output_file");

sub ofg_lineage {
    my $taxon = shift;
    my @lineage = $taxon->lineage;
    my %keep = map { $_ => 1 } qw[ order family genus ];
    return grep { $keep{$_->level} } reverse @lineage;
}
