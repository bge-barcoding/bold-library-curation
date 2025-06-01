#!/usr/bin/env perl
use strict;
use warnings;

# Script to fix column alignment issues in BAGS output
# Handles missing taxonomic levels by ensuring consistent column positions

my $input_file = shift @ARGV or die "Usage: $0 input_file [output_file]\n";
my $output_file = shift @ARGV || $input_file . ".fixed";

print STDERR "Fixing column alignment in: $input_file\n";
print STDERR "Output will be written to: $output_file\n";

open my $in_fh, '<', $input_file or die "Cannot open input file '$input_file': $!\n";
open my $out_fh, '>', $output_file or die "Cannot open output file '$output_file': $!\n";

my $line_count = 0;
my $fixed_count = 0;

while (my $line = <$in_fh>) {
    chomp $line;
    $line_count++;
    
    # Handle header line
    if ($line_count == 1) {
        print $out_fh $line . "\n";
        next;
    }
    
    my @fields = split /\t/, $line;
    my $field_count = scalar @fields;
    
    # Expected columns: taxonid, order, family, genus, species, BAGS, BIN, sharers
    if ($field_count == 8) {
        # Already correct, pass through
        print $out_fh $line . "\n";
    } elsif ($field_count == 7) {
        # Missing one taxonomic level, need to determine which one
        my ($taxonid, $bags_grade, $bin_url, $sharers);
        my @taxonomy;
        
        # Extract known fields
        $taxonid = $fields[0];
        $bags_grade = $fields[-3];  # BAGS grade is 3rd from end
        $bin_url = $fields[-2];     # BIN URL is 2nd from end  
        $sharers = $fields[-1];     # Sharers is last
        
        # Extract taxonomy fields (everything between taxonid and BAGS grade)
        @taxonomy = @fields[1..($field_count-4)];
        
        # Determine which level is missing and fix it
        my ($order, $family, $genus, $species);
        
        if (scalar @taxonomy == 3) {
            # One taxonomic level missing - need to determine which
            # Check if last taxonomy field looks like species name (has spaces or lowercase)
            if ($taxonomy[-1] =~ /\s/ || $taxonomy[-1] =~ /[a-z]/) {
                # Last field is species, so we have family, genus, species (missing order)
                $order = '';
                $family = $taxonomy[0];
                $genus = $taxonomy[1];
                $species = $taxonomy[2];
            } else {
                # Assume we have order, family, genus (missing species)
                # This shouldn't happen in species-level analysis, but handle it
                $order = $taxonomy[0];
                $family = $taxonomy[1];
                $genus = $taxonomy[2];
                $species = '';
            }
        } elsif (scalar @taxonomy == 2) {
            # Two taxonomic levels missing
            # Assume we have family, genus, species (missing order)
            $order = '';
            $family = $taxonomy[0];
            $genus = $taxonomy[1];
            $species = '';
        } else {
            # Unexpected case - just pad with empty strings
            $order = '';
            $family = $taxonomy[0] || '';
            $genus = $taxonomy[1] || '';
            $species = $taxonomy[2] || '';
        }
        
        # Reconstruct the line with proper column alignment
        my @fixed_fields = ($taxonid, $order, $family, $genus, $species, $bags_grade, $bin_url, $sharers);
        print $out_fh join("\t", @fixed_fields) . "\n";
        $fixed_count++;
        
    } elsif ($field_count == 6) {
        # Missing two taxonomic levels
        my ($taxonid, $bags_grade, $bin_url, $sharers);
        my @taxonomy;
        
        $taxonid = $fields[0];
        $bags_grade = $fields[-3];
        $bin_url = $fields[-2];
        $sharers = $fields[-1];
        
        @taxonomy = @fields[1..($field_count-4)];
        
        # Assume we have genus, species (missing order and family)
        my $order = '';
        my $family = '';
        my $genus = $taxonomy[0] || '';
        my $species = $taxonomy[1] || '';
        
        my @fixed_fields = ($taxonid, $order, $family, $genus, $species, $bags_grade, $bin_url, $sharers);
        print $out_fh join("\t", @fixed_fields) . "\n";
        $fixed_count++;
        
    } else {
        # Unexpected number of fields - print warning and pass through
        print STDERR "Warning: Line $line_count has unexpected field count ($field_count): $line\n";
        print $out_fh $line . "\n";
    }
}

close $in_fh;
close $out_fh;

print STDERR "Processing complete:\n";
print STDERR "  Total lines processed: $line_count\n";
print STDERR "  Lines fixed: $fixed_count\n";
print STDERR "  Output file: $output_file\n";

if ($fixed_count > 0) {
    print STDERR "\nFixed file is ready. You may want to:\n";
    print STDERR "1. Backup the original file\n";
    print STDERR "2. Replace the original with the fixed version\n";
    print STDERR "3. Re-run the import step\n";
}
