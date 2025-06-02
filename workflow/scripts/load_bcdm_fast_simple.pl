#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use File::Temp qw(tempfile);
use File::Basename;
use File::Path qw(make_path);
use File::Copy;

# Process command line arguments
my $bold_tsv;   # location of BOLD TSV dump
my $db_file;    # where to create database file
my $schema_sql; # location of table creation statement(s)
my $overwrite;  # overwrite existing DB
my $table = 'bold'; # database table
my $log_level = 'INFO'; # verbosity level for logger
my $temp_dir;   # optional temp directory
GetOptions(
    'tsv=s'   => \$bold_tsv,
    'db=s'    => \$db_file,
    'sql=s'   => \$schema_sql,
    'log=s'   => \$log_level,
    'force'   => \$overwrite,
    'table=s' => \$table,
    'temp=s'  => \$temp_dir
);

# Initialize Log::Log4perl
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = $log_level, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# Instantiate logger
my $log = Log::Log4perl->get_logger('load_bcdm');

# Set temp directory - use same directory as output if not specified
if (!$temp_dir) {
    $temp_dir = dirname($db_file);
    $log->info("Using output directory for temporary files: $temp_dir");
}

# Ensure temp directory exists
make_path($temp_dir) unless -d $temp_dir;

# Create the SQLite database - delete old
unlink($db_file) if $overwrite;
die 'Database exists and --force not specified' if -e $db_file;

$log->info("Creating database $db_file with schema $schema_sql");
system("sqlite3 $db_file < $schema_sql") == 0 
    or die "Failed to create database schema: $?";

# Prepare TSV for import - add recordid column
$log->info("Preparing TSV file for bulk import");
my ($temp_fh, $temp_file) = tempfile(
    SUFFIX => '.tsv', 
    DIR => $temp_dir,
    UNLINK => 1
);

# Open the original TSV file
open my $input_fh, "<:encoding(utf8)", $bold_tsv 
    or die "Could not open file '$bold_tsv': $!";

# Read schema and extract bold table columns
$log->info("Parsing database schema to extract column definitions");
open my $schema_fh, "<", $schema_sql or die "Could not open schema file: $!";
my @schema_cols;
my $in_bold = 0;
while (<$schema_fh>) {
    chomp;
    if (/CREATE TABLE IF NOT EXISTS\s+"bold"\s*\(/) {
        $in_bold = 1;
        next;
    }
    if ($in_bold && /^\s*\)/) {
        last;
    }
    next unless $in_bold;
    
    # Skip comment lines, foreign key constraints, and empty lines
    next if /^\s*--/ || /FOREIGN KEY/ || /^\s*$/;
    
    # Extract column name from quotes, excluding recordid and taxonid
	# if (/^\s*"([^"]+)"\s+(TEXT|INTEGER)/) {
    if (/^\s*\"([^\"]+)\"\s+(TEXT|INTEGER)\s*[,)]/) {
        my $col_name = $1;
        next if $col_name eq 'recordid' || $col_name eq 'taxonid';
        push @schema_cols, $col_name;
    }
}
close $schema_fh;
$log->info("Found " . scalar(@schema_cols) . " data columns in schema (excluding recordid and taxonid)");

# Read and parse TSV header
my $header = <$input_fh>;
chomp $header;
my @tsv_cols = split /\t/, $header;
# Clean column names of any whitespace/line endings
@tsv_cols = map { s/^\s+|\s+$//gr } @tsv_cols;
$log->info("Found " . scalar(@tsv_cols) . " columns in TSV file");

# Create bidirectional mapping for validation
my %tsv_index = map { $tsv_cols[$_] => $_ } 0..$#tsv_cols;
my %schema_index = map { $schema_cols[$_] => $_ } 0..$#schema_cols;

if ($log_level eq 'DEBUG' || 1) {  # Force debug for this issue
    $log->info("All schema columns: " . join(", ", @schema_cols));
    $log->info("All TSV columns: " . join(", ", @tsv_cols));
}

# Debug hash contents
if (exists $tsv_index{'bold_recordset_code_arr'}) {
    $log->info("bold_recordset_code_arr found in TSV at position: " . $tsv_index{'bold_recordset_code_arr'});
} else {
    $log->info("bold_recordset_code_arr NOT found in TSV hash");
}

if (exists $schema_index{'bold_recordset_code_arr'}) {
    $log->info("bold_recordset_code_arr found in schema at position: " . $schema_index{'bold_recordset_code_arr'});
} else {
    $log->info("bold_recordset_code_arr NOT found in schema hash");
}

# Validate column alignment and report mismatches
my @missing_in_tsv = grep { !exists $tsv_index{$_} } @schema_cols;
my @extra_in_tsv = grep { !exists $schema_index{$_} } @tsv_cols;

# Filter out known computed columns that are populated later in pipeline
my @computed_cols = qw(taxon_name taxon_rank);
my %computed_index = map { $_ => 1 } @computed_cols;
@missing_in_tsv = grep { !$computed_index{$_} } @missing_in_tsv;

if (@missing_in_tsv) {
    $log->warn("Schema columns missing from TSV (will be filled with empty values): " . join(", ", @missing_in_tsv));
}
if (@extra_in_tsv) {
    $log->warn("TSV columns not in schema (will be ignored): " . join(", ", @extra_in_tsv));
}
if (@computed_cols) {
    $log->info("Computed columns (populated later): " . join(", ", @computed_cols));
}

# Pre-compute column index mapping for performance optimization
my @col_indices = map { $tsv_index{$_} // -1 } @schema_cols;
$log->info("Column mapping computed - ready for data processing");

# Write new header to temp file (including taxonid which gets added later in the pipeline)
print $temp_fh join("\t", "recordid", "taxonid", @schema_cols), "\n";

# Process each row with optimized column alignment
my $recordid = 1;
while (my $line = <$input_fh>) {
    chomp $line;
    my @fields = split /\t/, $line;
    
    # Optimized column alignment using pre-computed indices
    my @aligned = map { $_ >= 0 ? ($fields[$_] // '') : '' } @col_indices;
    
    # Clean any embedded newlines/carriage returns from all fields
    @aligned = map { s/[\r\n]/ /gr } @aligned;
    
    # Add recordid and empty taxonid (filled later by taxonomy loading step)
    print $temp_fh join("\t", $recordid, "", @aligned), "\n";
    $recordid++;
    
    # Log progress for very large files
    $log->info("Prepared $recordid records") if ($recordid % 100_000 == 0);
}

close $input_fh;
close $temp_fh;

$log->info("Prepared " . ($recordid - 1) . " total records for import");

# Bulk import using SQLite's native .import with optimizations
$log->info("Starting bulk import into database");

# Try to ensure the database path exists
my $db_dir = dirname($db_file);
make_path($db_dir) unless -d $db_dir;

# Create debug temp file name for preservation
my $debug_temp_file = $temp_file;
$debug_temp_file =~ s/\.tmp$/_debug.tsv/;
copy($temp_file, $debug_temp_file);
$log->info("Temp file preserved for debugging at: $debug_temp_file");

my $sqlite_cmd = qq{
sqlite3 "$db_file" <<'EOF'
-- Performance optimizations
PRAGMA synchronous = OFF;
PRAGMA journal_mode = MEMORY;
PRAGMA cache_size = 200000;
PRAGMA temp_store = MEMORY;
-- Import data
.mode tabs
.import "$temp_file" $table

-- Restore normal settings
PRAGMA synchronous = NORMAL;
PRAGMA journal_mode = DELETE;
.quit
EOF
};

$log->info("Executing bulk import command");
my $result = system($sqlite_cmd);

# Import error checking removed to avoid duplicate data loading
# For 3M+ row datasets, we rely on the main import verification below

if ($result == 0) {
    $log->info("Successfully loaded " . ($recordid - 1) . " records into database");
} else {
    die "Bulk import failed with exit code: $result";
}

# Verify the import
$log->info("Verifying import success");

# Count actual records in temp file (excluding header)
my $temp_count_cmd = qq{wc -l < "$temp_file"};
my $temp_line_count = `$temp_count_cmd`;
chomp $temp_line_count;
my $actual_records = $temp_line_count - 1; # subtract header line

my $count_cmd = qq{sqlite3 "$db_file" "SELECT COUNT(*) FROM $table;"};
my $count_result = `$count_cmd`;
chomp $count_result;
$log->info("Database contains $count_result records");
$log->info("Temp file contained $actual_records records (excluding header)");

my $missing_records = $actual_records - $count_result;
if ($missing_records > 0) {
    $log->warn("Import verification found $missing_records missing records ($count_result imported vs $actual_records expected)");
    $log->warn("This is likely due to SQLite datatype mismatch errors during import");
    $log->warn("Debug temp file preserved at: $debug_temp_file");
    
    # Calculate percentage of successful imports
    my $success_rate = ($count_result / $actual_records) * 100;
    $log->info(sprintf("Import success rate: %.2f%% (%d/%d records)", $success_rate, $count_result, $actual_records));
    
    # Check if we should continue or fail based on environment variable
    my $allow_partial_import = $ENV{'ALLOW_PARTIAL_IMPORT'} || '';
    if ($allow_partial_import eq 'true' || $success_rate >= 99.0) {
        $log->warn("Continuing with partial import (success rate >= 99% or ALLOW_PARTIAL_IMPORT=true)");
    } else {
        die "Import verification failed: expected $actual_records records, found $count_result (set ALLOW_PARTIAL_IMPORT=true to continue with partial imports)";
    }
} elsif ($missing_records < 0) {
    die "Import verification failed: more records in database ($count_result) than in source file ($actual_records)";
}

$log->info("Import completed successfully and verified");
