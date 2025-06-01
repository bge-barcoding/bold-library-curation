use strict;
use warnings;
use Text::CSV;
use Getopt::Long;

# Cached BAGS rating, distinct sequences, file handles and taxonomic levels
my %BAGS;     # BIN=>BAGS/sharers mapping
my %SEEN;     # distinct sequences
my %HANDLE;   # file handles for families
my %SPECIES;  # species=>family mapping
my @LEVELS = qw[ phylum class order family genus species ]; # taxonomic levels

# Process command line arguments
my $dump_tsv; # where to access annotated BCDM dump file
my $bags_tsv; # location of BAGS TSV dump
my $fasta;    # where to write FASTA file
my $species;  # european species table, like in make_eu_specieslist.pl / european_pollinator_species.tsv
GetOptions(
    'dump=s'    => \$dump_tsv,
    'bags=s'    => \$bags_tsv,
    'fasta=s'   => \$fasta,
    'species=s' => \$species
);

# Initialize Text::CSV
my $tsv = Text::CSV->new({
    sep_char         => "\t",
    binary           => 1,
    auto_diag        => 1,
    allow_whitespace => 1,
    quote_char       => undef
});

# Create a hash table that maps BINs to BAGS ratings and sharers
{
    # Connect to the BAGS TSV as binary, with UTF-8, for reading
    open my $fh, "<:encoding(utf8)", $bags_tsv or die "Could not open file '$bags_tsv': $!";

    # Read the BAGS TSV as hash refs, store BIN=>BAGS mapping
    my @keys = @{$tsv->getline($fh)};
    $tsv->column_names(\@keys);
    while ( my $row = $tsv->getline_hr($fh) ) {

        # Cache BIN=>BAGS/sharers mapping - the BIN is the last part of the cluster URI in the TSV
        my $bin_url = $row->{'BIN'};
        if ( $bin_url =~ /clusteruri=(BOLD:.+)$/ ) {
            my $bin = $1;
            $BAGS{$bin} = {
                BAGS    => $row->{'BAGS'},
                sharers => $row->{'sharers'},
            };
        }
    }
}

# Create a hash table that maps species to families
{
    # Connect to the species TSV as binary, with UTF-8, for reading
    open my $fh, "<:encoding(utf8)", $species or die "Could not open file '$species': $!";

    # Read the species TSV as hash refs, store species=>family mapping
    my @keys = @{$tsv->getline($fh)};
    $tsv->column_names(\@keys);
    while ( my $row = $tsv->getline_hr($fh) ) {

        # Cache species=>family mapping
        $SPECIES{$row->{'species'}} = $row->{'family'};
    }
}

# Triage the BCDM dump TSV, writing to FASTA if BAGS=(A,B,C) and ranking=(1,2,3), otherwise to family TSVs
{
    # Connect to the BCDM TSV as binary, with UTF-8, for reading
    open my $fh, "<:encoding(utf8)", $dump_tsv or die "Could not open file '$dump_tsv': $!";

    # Open FASTA file for writing
    open my $fasta_fh, ">", $fasta  or die "Could not open file '$fasta': $!";

    # Read the BCDM dump TSV as hash refs
    my @keys = @{$tsv->getline($fh)};
    $tsv->column_names(\@keys);
    RECORD: while ( my $row = $tsv->getline_hr($fh) ) {

        # Skip if not having a species name or not in the species-to-family mapping hash
        next RECORD if not $row->{'species'} or not $SPECIES{$row->{'species'}};

        # Skip if not COI-5P
        next RECORD if not $row->{'marker_code'} or $row->{'marker_code'} ne 'COI-5P';

        # Skip if already seen the haplotype (sans gaps)
        my $seq = $row->{'nuc'};
        $seq =~ s/-//g;
        next RECORD if $SEEN{$seq}++;

        # Skip undefined or garbled BINs
        my $bin = $row->{'bin_uri'};
        next RECORD if not $bin or $bin !~ /^BOLD:.+$/;

        # BAGS grading > ABC, ranking > 3? Forward to family TSV with BAGS, record, ranking for manual curation
        if ( not $BAGS{$bin} or $BAGS{$bin}->{BAGS} !~ /^[ABC]$/ or not $row->{'ranking'} or $row->{'ranking'} > 3 ) {

            # Keep a pool of 1,000 handles. ulimit -Sn is 1024 on my system.
            if ( scalar(keys(%HANDLE)) == 1000 ) {
                my ($delete) = keys(%HANDLE);
                close $HANDLE{$delete};
                delete $HANDLE{$delete};
            }

            # Create or lookup file handle for family level TSV
            my $family = $SPECIES{$row->{'species'}};
            if ( not $HANDLE{$family} ) {

                # Print header upon file creation
                if ( not -e "$family.tsv" ) {
                    open my $family_fh, ">", "$family.tsv" or die "Could not open file '$family.tsv': $!";
                    print $family_fh join("\t", 'BAGS', 'sharers', @keys), "\n";
                    $HANDLE{$family} = $family_fh;
                }

                # Append to existing file
                else {
                    open my $family_fh, ">>", "$family.tsv" or die "Could not open file '$family.tsv': $!";
                    $HANDLE{$family} = $family_fh;
                }
            }
            my $family_fh = $HANDLE{$family};

            # Print the record with BAGS rating to family TSV
            no warnings 'uninitialized';
            print $family_fh join("\t", $BAGS{$bin}->{BAGS}, $BAGS{$bin}->{sharers}, map {$row->{$_}} @keys), "\n";
        }

        # Good BAGS rating
        else {

            # Print the record with lineage to FASTA
            my $process_id = $row->{'processid'};
            my $defline = join "|", 'private_BOLD', $process_id, $row->{'species'}, map {$row->{$_}} @LEVELS;
            print $fasta_fh ">$defline\n$seq\n";
        }
    }
}