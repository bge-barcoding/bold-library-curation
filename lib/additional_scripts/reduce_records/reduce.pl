#!/usr/bin/perl -w
use strict;
use warnings;

my %avail;
open(SPE, "all_specs_and_syn.csv");
while (<SPE>) {
    my $line = $_;
    chomp $line;
    my @array = split(";", $line);
    foreach (@array) {
        $avail{$_} = 1;
    }
}
close SPE;

unlink "reduced_BOLD_Public.04-Apr-2025.tsv";
open(OUT, ">>reduced_BOLD_Public.04-Apr-2025.tsv");

open(DAT, "<BOLD_Public.04-Apr-2025.tsv");

# --- Add header line ---
my $header = <DAT>;
chomp $header;
print OUT "$header\n";

# --- Process rest of file ---
while (<DAT>) {
    my $line = $_;
    chomp $line;
    my @array = split("\t", $line);
    if (defined $avail{$array[21]}) {
        print OUT "$line\n";
    }
}

close OUT;
close DAT;
