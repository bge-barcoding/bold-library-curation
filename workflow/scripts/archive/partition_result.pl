use strict;
use warnings;
use Getopt::Long;
use Text::CSV;
use IO::Uncompress::Gunzip qw(gunzip $GunzipError);
use IO::File;
use Data::Dumper;

# This script filters and partitions a BCDM TSV file. The --filter argument specifies a TSV file
# whose column headers must correspond with the BCDM. Every row in that input file specifies a 
# filter criterion which, if matched, allows the record to be retained. In the simplest case,
# this can be used to define a species list (with header `species`), which are retained in the
# output. The --column argument specifies how to partition the output. For example, if `marker_code`
# is given as a value, the output is partitioned into separate files, one for each marker. The
# --database argument specifies the location of the BCDM TSV.

# Command line options
my ($filter_file, $database_file, $column);
GetOptions(
    'filter=s'   => \$filter_file,   # resources/group_for_curation.tsv
    'database=s' => \$database_file, # results/BOLD TSV file
    'column=s'   => \$column,        # e.g. 'family'
) or die "Usage: $0 --filter <filter_file> --database <database_file> --column <column_name>\n";

# Check required arguments
die "Usage: $0 --filter <filter_file> --database <database_file> --column <column_name>\n" unless $filter_file && $database_file && $column;

# Read filter file using Text::CSV configured to operate on TSV, reading it as binary for UTF-8
open my $fh_filter, '<', $filter_file or die "Could not open '$filter_file': $!\n";
my $csv = Text::CSV->new({ sep_char => "\t", binary => 1 });
my $headers = $csv->getline($fh_filter);
$csv->column_names(@$headers);
my @filter_predicates;
while ( my $row = $csv->getline_hr($fh_filter) ) {
    push @filter_predicates, $row;
}
close $fh_filter;

# Open database file (handle gzipped file)
my $fh_database;
if ( $database_file =~ /\.gz$/ ) {
    $fh_database = IO::Uncompress::Gunzip->new($database_file) or die "gunzip failed: $GunzipError\n";
} 
else {
    open $fh_database, '<', $database_file or die "Could not open '$database_file': $!\n";
}

# Read database file header
my $csv_db = Text::CSV->new({ sep_char => "\t", binary => 1 });
my $db_headers = $csv_db->getline($fh_database);
$csv_db->column_names(@$db_headers);

# Check if the column exists in the database header
die "Column '$column' does not exist in the database file\n" unless grep { $_ eq $column } @$db_headers;

# Print database header
my $header_line = join("\t", @$db_headers);

# Function to check if a database row matches any filter predicate
sub matches_filter {
    my ( $db_row ) = @_;
    for my $filter ( @filter_predicates ) {
        my $match = 1;
        for my $key ( keys %$filter ) {
            if ( !defined $db_row->{$key} || $db_row->{$key} ne $filter->{$key} ) {
                $match = 0;
                last;
            }
        }
        return 1 if $match;
    }
    return 0;
}

# Hash to store file handles for each distinct column value
my %file_handles;

# Filter and print database rows
while ( my $db_row = $csv_db->getline_hr($fh_database) ) {
    if ( matches_filter($db_row) ) {
        my $column_value = $db_row->{$column};
        if ( !exists $file_handles{$column_value} ) {
            my $file_name = "${column_value}.tsv";
            my $fh = IO::File->new(">$file_name") or die "Could not open file '$file_name' for writing: $!\n";
            $file_handles{$column_value} = $fh;
            print $fh "$header_line\n";
        }
        my $fh = $file_handles{$column_value};
        print $fh join("\t", @{$db_row}{@$db_headers}), "\n";
    }
}

# Close all file handles
for my $fh ( values %file_handles ) {
    $fh->close;
}

close $fh_database;
