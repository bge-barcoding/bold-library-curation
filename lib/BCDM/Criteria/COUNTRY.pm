package BCDM::Criteria::COUNTRY;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::COUNTRY }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $country = $record->country_ocean;
    return ($country eq '' || $country eq 'None') ? 0 : 1, "Determined from country column";
}

1;
