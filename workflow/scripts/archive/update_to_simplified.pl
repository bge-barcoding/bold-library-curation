#!/usr/bin/env perl
use strict;
use warnings;

# Script to backup original assess_taxa_filtered.pl and replace with simplified version
# This allows immediate use of the simplified approach

my $script_dir = "C:\\GitHub\\bold-library-curation\\workflow\\scripts";
my $original_script = "$script_dir\\assess_taxa_filtered.pl";
my $simplified_script = "$script_dir\\assess_taxa_simplified.pl";
my $backup_script = "$script_dir\\assess_taxa_filtered.pl.backup";

print "Updating BAGS analysis script to use simplified output...\n\n";

# Check if files exist
unless (-f $original_script) {
    die "Error: Original script not found at $original_script\n";
}

unless (-f $simplified_script) {
    die "Error: Simplified script not found at $simplified_script\n";
}

# Create backup of original
if (-f $backup_script) {
    print "Backup already exists at $backup_script\n";
} else {
    print "Creating backup: $original_script -> $backup_script\n";
    system("copy", "\"$original_script\"", "\"$backup_script\"") == 0 
        or die "Failed to create backup: $!\n";
}

# Replace original with simplified version
print "Replacing original with simplified version...\n";
system("copy", "\"$simplified_script\"", "\"$original_script\"") == 0 
    or die "Failed to replace original script: $!\n";

print "\nUpdate completed successfully!\n";
print "- Original script backed up to: $backup_script\n";
print "- Simplified script is now active at: $original_script\n";
print "\nThe updated script outputs only 4 columns:\n";
print "  1. taxonid\n";
print "  2. BAGS (grade A-F)\n"; 
print "  3. BIN (URL)\n";
print "  4. sharers (comma-separated taxa)\n";
print "\nThis eliminates the column alignment issues caused by missing taxonomic ranks.\n";
