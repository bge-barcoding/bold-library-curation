package BCDM::Criteria::INSTITUTION;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::INSTITUTION }

my @neg = (
    'genbank',
    'no voucher',
    'none',
    'personal',
    'private',
    'research collection of',
    'unknown',
    'unvouchered'
);

my $log = __PACKAGE__->_get_logger(__PACKAGE__, 'DEBUG');

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
# Here, the criterion to assess is:
# 'Specimen is vouchered in an institution'
# This will involve substring matching in
# the inst column. There are many
# distinct values but there is some
# regularity to them.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $method = $record->inst;
    my $id = $record->recordid;

    # Check negative matches
    my @mn;
    for my $pattern ( @neg ) {
         if ( $method =~ /$pattern/i ) {
            push @mn, $pattern;
            $log->info("Negative match for $id: $pattern")
        }
    }

    # assume it passed
    my $result = 1;

    # failure if list of matches is non-empty
    if ( @mn ) {
        $result = 0;
    }

    # failure if field is empty
    elsif ( not $method ) {
        $result = 0;
    }

    # Return result
    return $result, "Based on whether it has a negative (@mn) or no matches (null)";
    
}

1;
