#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use File::Temp qw(tempfile);
use File::Copy;

# Process command line arguments
my $bold_tsv;   # location of BOLD TSV dump
my $db_file;    # where to create database file
my $schema_sql; # location of table creation statement(s)
my $overwrite;  # overwrite existing DB
my $table = 'bold'; # database table
my $log_level = 'INFO'; # verbosity level for logger
GetOptions(
    'tsv=s'   => \$bold_tsv,
    'db=s'    => \$db_file,
    'sql=s'   => \$schema_sql,
    'log=s'   => \$log_level,
    'force'   => \$overwrite,
    'table=s' => \$table
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

# Create the SQLite database - delete old
unlink($db_file) if $overwrite;
die 'Database exists and --force not specified' if -e $db_file;

$log->info("Creating database $db_file with schema $schema_sql");
system("sqlite3 $db_file < $schema_sql") == 0 
    or die "Failed to create database schema: $?";

# Prepare TSV for import - add recordid column
$log->info("Preparing TSV file for bulk import");
my ($temp_fh, $temp_file) = tempfile(SUFFIX => '.tsv', UNLINK => 1);

# Read original TSV and add recordid column
open my $input_fh, "<:encoding(utf8)", $bold_tsv 
    or die "Could not open file '$bold_tsv': $!";

# Read and modify header
my $header = <$input_fh>;
chomp $header;
print $temp_fh "recordid\t$header\n";

# Add recordid to each data row
my $recordid = 1;
while (my $line = <$input_fh>) {
    chomp $line;
    print $temp_fh "$recordid\t$line\n";
    $recordid++;
    
    # Log progress for very large files
    $log->info("Prepared $recordid records") if ($recordid % 100_000 == 0);
}

close $input_fh;
close $temp_fh;

$log->info("Prepared $recordid total records for import");

# Bulk import using SQLite's native .import with optimizations
$log->info("Starting bulk import into database");
my $sqlite_cmd = qq{
sqlite3 "$db_file" <<'EOF'
-- Performance optimizations
PRAGMA synchronous = OFF;
PRAGMA journal_mode = MEMORY;
PRAGMA cache_size = 200000;
PRAGMA temp_store = MEMORY;
PRAGMA mmap_size = 2147483648;

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

if ($result == 0) {
    $log->info("Successfully loaded $recordid records into database $db_file");
} else {
    die "Bulk import failed with exit code: $result";
}

# Verify the import
$log->info("Verifying import success");
my $count_result = `sqlite3 "$db_file" "SELECT COUNT(*) FROM $table;"`;
chomp $count_result;
$log->info("Database contains $count_result records");

if ($count_result != $recordid - 1) {  # -1 because recordid started at 1
    die "Import verification failed: expected " . ($recordid - 1) . " records, found $count_result";
}

$log->info("Import completed successfully and verified");
