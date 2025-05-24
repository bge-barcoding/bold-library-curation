package BCDM::Criteria::IDENTIFIER;
use strict;
use warnings;
use base 'BCDM::Criteria';

# values we assessed are not taxonomic experts (case-insensitive partial matches)
my @cbg_patterns = (
    'Kate Perez',
    'Angela Telfer',
    'BOLD',
    'BLAST',
    'BIN',
    'None',
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

    # check if identifier contains any of the non-expert patterns (case-insensitive)
    for my $pattern (@cbg_patterns) {
        if ( $identifier =~ /\Q$pattern\E/i ) {
            return 0, "identified_by: '$identifier' (contains '$pattern')";
        }
    }

    # everything else - assumed to be a taxonomic expert
    return 1, "identified_by: '$identifier'";
}

1;
