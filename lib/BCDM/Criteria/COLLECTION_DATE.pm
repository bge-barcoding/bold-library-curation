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
    
    my $start_date = $record->collection_date_start;
    my $end_date = $record->collection_date_end;
    
    # Check if we have valid dates (not undef, empty, or "None")
    if ( ($start_date && $start_date ne "None") or ($end_date && $end_date ne "None") ) {
        $has_date = 1;
    }
    return $has_date, "Determined from collection_date column";
}

1;
