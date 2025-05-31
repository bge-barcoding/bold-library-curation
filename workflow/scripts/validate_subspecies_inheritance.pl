#!/usr/bin/env perl
# Quick validation script for subspecies BAGS inheritance
use strict;
use warnings;
use DBI;

my $db_file = $ARGV[0] or die "Usage: $0 <database_file>\n";

my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", "", "", {
    RaiseError => 1,
    AutoCommit => 1
});

print "=== Subspecies BAGS Inheritance Validation ===\n\n";

# Count total subspecies
my $total_subspecies = $dbh->selectrow_array(
    "SELECT COUNT(*) FROM taxa WHERE level = 'subspecies'"
);
print "Total subspecies in database: $total_subspecies\n";

# Count subspecies with BAGS grades
my $subspecies_with_bags = $dbh->selectrow_array(
    "SELECT COUNT(DISTINCT s.taxonid) FROM taxa s 
     JOIN bags b ON s.taxonid = b.taxonid 
     WHERE s.level = 'subspecies'"
);
print "Subspecies with BAGS grades: $subspecies_with_bags\n";

if ($total_subspecies > 0) {
    my $percentage = sprintf("%.1f", ($subspecies_with_bags / $total_subspecies) * 100);
    print "Inheritance coverage: $percentage%\n\n";
}

# Show grade distribution for subspecies
print "Grade distribution for subspecies:\n";
my $grade_dist = $dbh->selectall_arrayref(
    "SELECT b.bags_grade, COUNT(*) as count
     FROM taxa s JOIN bags b ON s.taxonid = b.taxonid
     WHERE s.level = 'subspecies'
     GROUP BY b.bags_grade ORDER BY b.bags_grade"
);

foreach my $row (@$grade_dist) {
    printf "  Grade %s: %d subspecies\n", $row->[0], $row->[1];
}

print "\nValidation complete!\n";
$dbh->disconnect;
