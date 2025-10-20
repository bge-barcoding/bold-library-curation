package BCDM::Criteria::COLLECTORS;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::COLLECTORS }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $collectors = $record->collectors;
    return ($collectors eq '' || lc($collectors) eq 'none') ? 0 : 1, "Determined from collectors column";
}

1;
