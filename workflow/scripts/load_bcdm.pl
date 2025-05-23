use strict;
use warnings;
use DBI;
use Text::CSV;
use Getopt::Long;
use Log::Log4perl qw(:easy);

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
$log->info("Going to create database $db_file with schema $schema_sql");
system("sqlite3 $db_file < $schema_sql");

# Connect to it with the DataBase Interface (DBI) module system
$log->info("Going to connect to created database $db_file using the DBI");
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file","","");

# Connect to the BOLD TSV as binary, with UTF-8
$log->info("Going to read BOLD TSV dump $bold_tsv");
my $tsv = Text::CSV->new({ 
    sep_char => "\t", 
    binary => 1, 
    auto_diag => 1, 
    allow_whitespace => 1, 
    quote_char => '"',
    escape_char => '"',
    allow_loose_quotes => 1,
    allow_loose_escapes => 1
});
open my $fh, "<:encoding(utf8)", $bold_tsv or die "Could not open file '$bold_tsv': $!";

# Prepare the SQL statement
my @keys = @{ $tsv->getline($fh) };
my $columns = join ", ", 'recordid', map { $dbh->quote($_) } @keys;
my $count = scalar(@keys) + 1;
my $qmarks = join ", ", ('?') x $count;
my $sql = "INSERT INTO $table ($columns) VALUES ($qmarks)";
my $sth = $dbh->prepare($sql) or die $dbh->errstr;

# Iterate over the rows as hash refs
my $PK = 1;
$tsv->column_names(\@keys);
while (my $row = $tsv->getline_hr($fh)) {

    # Insert row values prefixed with primary key
    my @values = ( $PK++, map { $row->{$_} } @keys );
    $sth->execute(@values) or die $dbh->errstr;

    # Log progress
    $log->info("Processed record $PK") unless $PK % 10_000;
}
