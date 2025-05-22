package BCDM::Criteria;
use strict;
use warnings;
use Carp 'croak';
use BCDM::ORM;
use Module::Load;
use Log::Log4perl qw(:easy);

our $SPECIES_ID=1;      # 455.63s user 8.64s system 95% cpu 8:05.72 total
our $TYPE_SPECIMEN=2;   # 1411.82s user 20.51s system 93% cpu 25:24.25 total
our $SEQ_QUALITY=3;     # 846.79s user 10.82s system 97% cpu 14:37.19 total
our $PUBLIC_VOUCHER=4;  # 1084.40s user 17.85s system 96% cpu 19:01.70 total
our $HAS_IMAGE=5;       # 2553.98s user 208.97s system 3% cpu 25:06:17.02 total
our $IDENTIFIER=6;      # 440.84s user 7.14s system 97% cpu 7:41.26 total
our $ID_METHOD=7;       # 874.82s user 20.55s system 94% cpu 15:44.43 total
our $COLLECTORS=8;      # 484.61s user 11.80s system 95% cpu 8:40.16 total
our $COLLECTION_DATE=9; # 516.87s user 14.21s system 93% cpu 9:27.06 total
our $COUNTRY=10;        # 454.66s user 8.69s system 96% cpu 8:01.04 total
our $SITE=11;           # 466.81s user 9.39s system 94% cpu 8:21.46 total
our $REGION=12;
our $SECTOR=13;
our $COORD=14;          # 451.98s user 7.73s system 96% cpu 7:55.33 total
our $INSTITUTION=15;
our $MUSEUM_ID=16;

# Initialize Log::Log4perl
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = INFO, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# input is a bold table record and a
# string version of the name, e.g. 'SPECIES_ID'
{
    my @queue;
    sub assess {
        my ($self, %args) = @_;

        # add the record to the queue
        push @queue, $args{record};

        # delegate the assessment if the queue is full
        if ( scalar(@queue) == $self->_batch_size ) {
            my @result = $self->_assess(@queue);
            for ( my $i = 0; $i <= $#result - 1; $i += 2 ) {
                my $status = $result[$i];
                my $notes  = $result[$i+1];
                $args{handler}->( $status, $notes, (($i/2)+1) );
            }
            @queue = ();
        }
    }
}

sub persist {
    my ( $self, %args ) = @_;
    my $cid    = $self->_criterion;
    my $record = $args{'record'};
    my $status = $args{'status'};
    my $notes  = $args{'notes'};
    my $schema = $record->result_source->schema;
    my $result = $schema->resultset('BoldCriteria')->find_or_create({
        criterionid => $cid,
        recordid    => $record->recordid,
    });
    $result->update({ status => $status, notes => $notes });
}

sub _get_logger {
    my ( $self, $name ) = @_;
    return Log::Log4perl->get_logger($name);
}

sub load_criterion {
    my ( $self, $criterion ) = @_;
    my $package = __PACKAGE__ . '::' . uc($criterion);
    eval { load $package };
    if ( $@ ) {
        croak "Unknown criterion $criterion: $@";
    }
    return $package;
}

sub _batch_size { 1 }

1;
