package BCDM::Criteria::COORD;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::COORD }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $coord = $record->coord;
    return ($coord eq '' || lc($coord) eq 'none') ? 0 : 1, "Determined from coord column";
}

1;
