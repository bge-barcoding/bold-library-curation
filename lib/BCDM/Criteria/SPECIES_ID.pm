package BCDM::Criteria::SPECIES_ID;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::SPECIES_ID }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
# Here, the criterion to assess is:
# 'Specimen was identified to species rank'
sub _assess {
    my $package = shift;
    my $record = shift;
    
# Check if the species field contains excluded patterns, is empty, or is "None"
    if ($record->species =~ /^\s*$/ || 
        $record->species =~ /\s+cf\s+/ || 
        $record->species =~ /\s+nr\s+/ || 
        $record->species =~ /\s+aff\s+/ || 
        $record->species =~ /\d+/ || 
        $record->species =~ /\./ || 
        $record->species =~ /^\s*None\s*$/i) {
        return 0, "Species column is empty, contains excluded patterns ( cf , nr , aff , numbers, or periods), or is 'None'";
    }

# If species field is not empty and doesn't contain excluded patterns or 'None', the criterion passes
    return 1, "Determined from species column";
}
    
1;
