#!/usr/bin/env perl
use strict;
use warnings;

# Script to convert existing BAGS output to simplified format
# Extracts only taxonid, BAGS, BIN, and sharers columns

my $input_file = shift @ARGV or die "Usage: $0 input_file [output_file]\n";
my $output_file = shift @ARGV || $input_file . ".simplified";

print STDERR "Simplifying BAGS output: $input_file\n";
print STDERR "Output will be written to: $output_file\n";

open my $in_fh, '<', $input_file or die "Cannot open input file '$input_file': $!\n";
open my $out_fh, '>', $output_file or die "Cannot open output file '$output_file': $!\n";

my $line_count = 0;
my $processed_count = 0;

# Output simplified header
print $out_fh join("\t", qw[taxonid BAGS BIN sharers]), "\n";

while (my $line = <$in_fh>) {
    chomp $line;
    $line_count++;
    
    # Skip original header line
    if ($line_count == 1) {
        next;
    }
    
    my @fields = split /\t/, $line;
    my $field_count = scalar @fields;
    
    # Extract the essential columns regardless of current alignment issues
    my $taxonid = $fields[0];
    my $bags_grade = '';
    my $bin_url = '';
    my $sharers = '';
    
    # Find BAGS grade, BIN URL, and sharers from the end of the line
    if ($field_count >= 3) {
        # BAGS grade is typically 3rd from end, BIN is 2nd from end, sharers is last
        # But we need to handle cases where columns are misaligned
        
        # Look for the BIN URL (contains 'boldsystems.org')
        my $bin_index = -1;
        for my $i (0..$#fields) {
            if ($fields[$i] =~ /boldsystems\.org/) {
                $bin_index = $i;
                last;
            }
        }
        
        if ($bin_index >= 0) {
            # Found BIN URL
            $bin_url = $fields[$bin_index];
            $sharers = $fields[$bin_index + 1] || '';
            
            # BAGS grade should be the field before BIN URL
            if ($bin_index > 0) {
                $bags_grade = $fields[$bin_index - 1];
            }
        } else {
            # Fallback: assume last 3 fields are BAGS, BIN, sharers
            $bags_grade = $fields[-3] || '';
            $bin_url = $fields[-2] || '';
            $sharers = $fields[-1] || '';
        }
    }
    
    # Validate that we have the essential data
    if ($taxonid && $bags_grade =~ /^[A-F]$/ && $bin_url =~ /boldsystems\.org/) {
        print $out_fh join("\t", $taxonid, $bags_grade, $bin_url, $sharers), "\n";
        $processed_count++;
    } else {
        print STDERR "Warning: Line $line_count could not be processed properly: $line\n" if $line_count <= 10;
    }
}

close $in_fh;
close $out_fh;

print STDERR "Processing complete:\n";
print STDERR "  Total lines processed: $line_count\n";
print STDERR "  Valid records output: $processed_count\n";
print STDERR "  Output file: $output_file\n";

print STDERR "\nSimplified file is ready with 4 columns:\n";
print STDERR "  1. taxonid\n";
print STDERR "  2. BAGS (grade)\n";
print STDERR "  3. BIN (URL)\n";
print STDERR "  4. sharers (comma-separated)\n";
