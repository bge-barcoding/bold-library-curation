package BCDM::Criteria::SEQ_QUALITY;
use strict;
use warnings;
use base 'BCDM::Criteria';

# this so that we know the criterionid for
# updates in the intersection table
sub _criterion { $BCDM::Criteria::SEQ_QUALITY }

# this tests the criterion and returns
# boolean 0/1 depending on fail/pass. In
# addition, optional notes may be returned.
# Here, the criterion to assess is:
# 'Sequence is long, with few ambiguities'
# This will involve aligning each sequence
# against an HMM (such as the one from
# FinPROTAX), trimming the ends outside the
# 658 bp range, then trimming the gap chars
# and ambiguities, and reporting the 
# remainder (notes?) and summarizing that
# as pass/fail.
sub _assess {
    my $package = shift;
    my $record  = shift;
    my $logger  = $package->_get_logger(__PACKAGE__, 'INFO');
    my $bin_uri = $record->bin_uri;
    my $nucraw  = $record->nuc;

    # remove leading and trailing non-ACGT characters
    $logger->debug('Removing leading/trailing non-ACGT characters');
    $nucraw =~ s/^[^ACGTacgt]*//;
    $nucraw =~ s/[^ACGTacgt]*$//;

    # count the number of unambiguous characters here
    my $count = ($nucraw =~ tr/ACGTacgt//);

    # TODO have minlength managed by config.yml
    if ( $bin_uri ne 'none' and $count >= 500 ) {
        return 1, 'Has BIN assignment and length >= 500bp';
    }
    else {
        return 0, 'No BIN assignment and/or length < 500bp';
    }
}

1;
