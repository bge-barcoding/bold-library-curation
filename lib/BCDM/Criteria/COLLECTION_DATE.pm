package BCDM::Criteria::COLLECTION_DATE;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::COLLECTION_DATE }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $has_date = 0;
    if ( $record->collection_date_start or $record->collection_date_end ) {
        $has_date = 1;
    }
    return $has_date, "Determined from collection_date column";
}

1;
