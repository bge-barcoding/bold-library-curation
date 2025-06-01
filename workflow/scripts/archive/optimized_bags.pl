package BCDM::BAGS::Optimized;
use strict;
use warnings;
use Carp;
use Data::Dumper;
use BCDM::Criteria;

our $AUTOLOAD;
my $log = BCDM::Criteria->get_logger('BAGS::Optimized');

# Class-level cache for BIN sharing data to avoid repeated queries
our %BIN_SHARING_CACHE;
our %SUBSPECIES_CACHE;

sub new {
    my ( $class, $taxon, $orm_handle ) = @_;
    my $self = {
        'taxon'     => undef,
        'bins'      => [],
        'records'   => [],
        'orm'       => $orm_handle,  # Keep ORM handle for reuse
    };
    bless $self, $class;
    $self->taxon($taxon);
    return $self;
}

sub taxon {
    my ( $self, $taxon ) = @_;
    if ( defined $taxon ) {
        my ( $id, $name, $level, $kingdom ) = map { $taxon->$_ } qw( taxonid name level kingdom );
        $log->debug("Setting $level $name ($id) from kingdom $kingdom");
        $self->{'taxon'} = $taxon;

        # Use prepared statement for better performance
        my $orm = $self->{'orm'} || $taxon->result_source->schema;
        
        # Optimize subspecies lookup with caching
        my @search_conditions = ({ taxonid => $id });
        
        if ( $level eq 'species' ) {
            my $subspecies_taxonids = $self->_get_subspecies_taxonids($orm, $name);
            if ( @$subspecies_taxonids ) {
                $log->debug("Including " . scalar(@$subspecies_taxonids) . " subspecies records for $name");
                push @search_conditions, { taxonid => { -in => $subspecies_taxonids } };
            }
        }
        
        # Use optimized query with proper indexing
        my $search_condition = @search_conditions > 1 ? 
            { -or => \@search_conditions } : 
            $search_conditions[0];
        
        # Use iterator to avoid loading all records into memory at once
        my $records_rs = $orm->resultset('Bold')->search( 
            $search_condition,
            { 
                columns => [qw/recordid bin_uri taxonid/],  # Only load needed columns
                result_class => 'DBIx::Class::ResultClass::HashRefInflator'
            }
        );
        
        my @records;
        my %bins_seen;
        while (my $record = $records_rs->next) {
            push @records, $record;
            $bins_seen{$record->{bin_uri}}++ if defined $record->{bin_uri};
        }
        
        $self->records( \@records );
        $self->bins( [ keys %bins_seen ] );
        
        $log->debug("Found " . $self->n_records . " records, " . $self->n_bins . " BINs for $name");
    }
    return $self->{'taxon'};
}

sub _get_subspecies_taxonids {
    my ($self, $orm, $species_name) = @_;
    
    # Use cache to avoid repeated subspecies lookups
    return $SUBSPECIES_CACHE{$species_name} if exists $SUBSPECIES_CACHE{$species_name};
    
    my $species_name_pattern = $species_name . ' %';
    my $subspecies_rs = $orm->resultset('Taxa')->search(
        {
            level => 'subspecies',
            name => { -like => $species_name_pattern }
        },
        {
            columns => ['taxonid'],
            result_class => 'DBIx::Class::ResultClass::HashRefInflator'
        }
    );
    
    my @subspecies_taxonids;
    while (my $subspecies = $subspecies_rs->next) {
        push @subspecies_taxonids, $subspecies->{taxonid};
    }
    
    $SUBSPECIES_CACHE{$species_name} = \@subspecies_taxonids;
    return \@subspecies_taxonids;
}

sub n_records { scalar @{ shift->records } }

sub n_bins { scalar @{ shift->bins } }

sub grade {
    my $self = shift;
    
    # Use cached sharing data for faster grade calculation
    my $sharing_species_count = scalar($self->sharing_taxa);
    my $is_shared = $sharing_species_count > 0;

    # Grade logic (same as original)
    if ( $self->n_records >= 10 && $self->n_bins == 1 && !$is_shared ) {
        return 'A';
    }
    elsif ( 3 <= $self->n_records && $self->n_records < 10 && $self->n_bins == 1 && !$is_shared ) {
        return 'B';
    }
    elsif ( $self->n_bins > 1 && !$is_shared ) {
        return 'C';
    }
    elsif ( $self->n_records < 3 && $self->n_bins == 1 && !$is_shared ) {
        return 'D';
    }
    elsif ($is_shared) {
        return 'E';
    }
    else {
        $log->warn("Could not determine BAGS grade for " . $self->taxon->name);
        return 'F';
    }
}

sub taxa_sharing_bin {
    my ( $self, $bin ) = @_;
    
    # Use cache to avoid repeated BIN sharing queries
    return @{$BIN_SHARING_CACHE{$bin}} if exists $BIN_SHARING_CACHE{$bin};
    
    my $orm = $self->{'orm'} || $self->taxon->result_source->schema;

    # Optimized query using joins instead of nested loops
    my $sharing_rs = $orm->resultset('Bold')->search(
        { bin_uri => $bin },
        {
            join => 'taxon',
            columns => ['taxon.name', 'taxon.level'],
            distinct => 1,
            result_class => 'DBIx::Class::ResultClass::HashRefInflator'
        }
    );
    
    my %seen_species;
    my $current_species = $self->taxon->name;
    
    while (my $record = $sharing_rs->next) {
        next unless defined $record->{name};
        
        my $level = $record->{level} // '';
        my $name = $record->{name};
        
        if ( $level eq 'species' ) {
            $seen_species{$name}++ if $name ne $current_species;
        } elsif ( $level eq 'subspecies' ) {
            # Extract parent species name
            my @name_parts = split /\s+/, $name;
            if (@name_parts >= 2) {
                my $parent_species_name = join(' ', @name_parts[0,1]);
                $seen_species{$parent_species_name}++ if $parent_species_name ne $current_species;
            }
        }
    }
    
    my @sharing_species = keys %seen_species;
    $BIN_SHARING_CACHE{$bin} = \@sharing_species;
    
    return @sharing_species;
}

sub sharing_taxa {
    my $self = shift;
    
    my %all_sharing_species;
    for my $bin ( @{ $self->bins } ) {
        next unless defined $bin && $bin =~ /^BOLD:/;
        my @sharing = $self->taxa_sharing_bin($bin);
        for my $species (@sharing) {
            $all_sharing_species{$species}++;
        }
    }
    
    return keys %all_sharing_species;
}

# Clear cache periodically to prevent memory buildup
sub clear_cache {
    %BIN_SHARING_CACHE = ();
    %SUBSPECIES_CACHE = ();
}

sub AUTOLOAD {
    my ( $self, $arg ) = @_;
    my $method = $AUTOLOAD;
    $method =~ s/.*:://;

    if ( exists $self->{$method} ) {
        if ( defined $arg ) {
            $self->{$method} = $arg;
        }
        return $self->{$method};
    }
    elsif ( $method =~ /^[A-Z]+$/ ) {
        return;
    }
    else {
        Carp::carp Dumper($self);
        Carp::croak "No such method '$method'";
    }
}

1;
