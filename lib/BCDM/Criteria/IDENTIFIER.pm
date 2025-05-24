package BCDM::Criteria::IDENTIFIER;
use strict;
use warnings;
use base 'BCDM::Criteria';

# values we assessed are not taxonomic experts
# Using word boundaries (\b) to avoid matching parts of legitimate names
my @cbg_patterns = (
    qr/\bKate\s+Perez\b/i,      # Exact name match - CBG employee
    qr/\bAngela\s+Telfer\b/i,   # Exact name match - CBG employee  
    qr/\bBOLD\b/i,              # Word boundary - won't match names containing "bold" like "Theobold"
    qr/\bBLAST\b/i,             # Word boundary - won't match names containing "blast" like "Blaster" ... names these days!
    qr/\bBIN\b/i,               # Word boundary - won't match names containing "bin" like "Sabine"
    qr/\bNone\b/i,              # Word boundary - won't match names containing "none" like "Nonel"
);

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::IDENTIFIER }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
# Here, the criterion to assess is:
# 'Specimen was identified by a named person'
sub _assess {
    my $package = shift;
    my $record = shift;
    my $identifier = $record->identified_by;

    # no identifier given
    if ( not $identifier ) {
        return 0, "no identifier named";
    }

    # check if identifier matches any of the non-expert patterns
    for my $pattern (@cbg_patterns) {
        if ( $identifier =~ /$pattern/ ) {
            return 0, "identified_by: '$identifier' (matches exclusion pattern)";
        }
    }

    # everything else - assumed to be a taxonomic expert
    return 1, "identified_by: '$identifier'";
}

1;
