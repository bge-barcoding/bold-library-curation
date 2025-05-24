package BCDM::Criteria::ID_METHOD;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::ID_METHOD }

my @pos = qw(
    descr
    det
    diss
    exam
    expert
    genit
    identifier
    key
    label
    literature
    micros
    mor
    taxonomic
    type
    vou
    guide
    flora
    specimen
    traditional
    visual
    wing
    logical
    knowledge
    photo
    verified
    key
);

my @neg = qw(
    barco
    BOLD
    mBRAVE
    SINTAX
    CO1
    COI
    COX
    DNA
    mole
    phylo
    sequ
    tree
    image
    bin
    silva
    ncbi
    ncbl
    engine
    blast
    genbank
    genetic
    unspecified
    its
);

my $log = __PACKAGE__->_get_logger(__PACKAGE__, 'DEBUG');

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
# Here, the criterion to assess is:
# 'Specimen was identified by morphology'
# This will involve substring matching against
# identification_method. It looks like there
# are many that start with '^BIN' or with
# '^BOLD'. Those are certainly reverse
# taxonomy.
sub _assess {
    my $package = shift;
    my $record = shift;
    my $method = $record->identification_method;
    my $id = $record->recordid;


    # Initialize flags for pattern matching
    my $has_positive = 0;
    my $has_negative = 0;
    
    # Check for positive matches (morphological identification indicators)
    for my $pattern ( @pos ) {
         if ( $method =~ /$pattern/i ) {
            $has_positive = 1;
            $log->info("Positive match for $id: $pattern");
            last; # Found at least one positive, no need to continue
        }
    }
    
    # Check for negative matches (molecular/automated identification indicators)
    for my $pattern ( @neg ) {
         if ( $method =~ /$pattern/i ) {
            $has_negative = 1;
            $log->info("Negative match for $id: $pattern");
            last; # Found at least one negative, no need to continue
        }
    }
    
    # Determine score based on pattern combinations:
    # 0 = No positive matches (includes no matches at all OR negative only)
    # 1 = Positive matches only (pure morphological identification)
    # 2 = Both positive and negative matches (mixed methods)
    my $retval;
    my $notes;
    
    if (!$has_positive) {
        # No positive matches found - either no matches or negative only
        $retval = 0;
        $notes = "No morphological identification indicators found";
    } elsif ($has_positive && !$has_negative) {
        # Positive matches only - pure morphological identification
        $retval = 1;
        $notes = "Pure morphological identification indicators found";
    } else {
        # Both positive and negative matches - mixed methods
        $retval = 2;
        $notes = "Mixed identification methods detected (both morphological and molecular indicators)";
    }
    
    $log->info("Final score for $id: $retval ($notes)");
    
    # Return result
    return $retval, $notes;
}

1;
