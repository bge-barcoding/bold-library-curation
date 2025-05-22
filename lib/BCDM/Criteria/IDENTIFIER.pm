package BCDM::Criteria::IDENTIFIER;
use strict;
use warnings;
use base 'BCDM::Criteria';

# values we assessed are not taxonomic experts
my %cbg = (
    'Kate Perez'     => 1,
    'Angela Telfer'  => 1,
    'BOLD' => 1,
    'BLAST' => 1,
    'BIN' => 1,
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

    # known entity not considered an expert
    if ( $cbg{$identifier} ) {
        return 0, "identified_by: '$identifier'";
    }

    # no identifier given
    elsif ( not $identifier ) {
        return 0, "no identifier named";
    }

    # everything else
    else {
        return 1, "identified_by: '$identifier'";
    }
}

1;
