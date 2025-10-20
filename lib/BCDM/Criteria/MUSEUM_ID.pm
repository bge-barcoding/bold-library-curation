package BCDM::Criteria::MUSEUM_ID;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::MUSEUM_ID }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $museumid = $record->museumid;
    return ($museumid eq '' || lc($museumid) eq 'none') ? 0 : 1, "Determined from museumid column";
}

1;
