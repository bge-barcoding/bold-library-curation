package BCDM::BAGS;
use strict;
use warnings;
use Carp;
use Data::Dumper;
use BCDM::Criteria;

our $AUTOLOAD;
my $log = BCDM::Criteria->get_logger('BAGS');

sub new {
    my ( $class, $taxon ) = @_;
    my $self = {
        'taxon'   => undef, # BCDM::ORM::Result::Taxa
        'bins'    => [],
        'records' => [],
    };
    bless $self, $class;
    $self->taxon($taxon);
    return $self;
}

sub taxon {
    my ( $self, $taxon ) = @_;
    if ( defined $taxon ) {
        my ( $id, $name, $level, $kingdom ) = map { $taxon->$_ } qw( taxonid name level kingdom );
        $log->info("Setting $level $name ($id) from kingdom $kingdom");
        $self->{'taxon'} = $taxon;

        # Get all barcodes for this taxon - REMOVED COI-5P filter as it's handled upstream
        # Include records that are identified to species level (same taxonid) 
        # AND subspecies that belong to this species
        my $orm = $taxon->result_source->schema;
        
        # Build the query to include both species-level and subspecies records
        my @search_conditions;
        
        # Always include records with the exact taxonid (species-level records)
        push @search_conditions, { taxonid => $id };
        
        # If this is a species-level taxon, also include subspecies records
        if ( $level eq 'species' ) {
            # Find all subspecies that belong to this species by looking for taxa where:
            # 1. Have level = 'subspecies' 
            # 2. The name starts with this species name (e.g., "Argia fumipennis violacea")
            my $species_name_pattern = $name . ' %';
            my $subspecies_rs = $orm->resultset('Taxa')->search({
                level => 'subspecies',
                name => { -like => $species_name_pattern }
            });
            
            my @subspecies_taxonids;
            while ( my $subspecies = $subspecies_rs->next ) {
                push @subspecies_taxonids, $subspecies->taxonid;
                $log->info("Found subspecies: " . $subspecies->name . " (taxonid " . $subspecies->taxonid . ")");
            }
            
            if ( @subspecies_taxonids ) {
                $log->info("Including " . scalar(@subspecies_taxonids) . " subspecies records for $name");
                push @search_conditions, { taxonid => { -in => \@subspecies_taxonids } };
            }
        }
        
        # Combine conditions with OR
        my $search_condition = @search_conditions > 1 ? 
            { -or => \@search_conditions } : 
            $search_conditions[0];
        
        $self->records( [ $orm->resultset('Bold')->search( $search_condition )->all ] );
        $log->info("Found " . $self->n_records . " species + subspecies records for $name");

        # Get all distinct, defined BINs for this taxon's records
        # Filter out undefined, empty, and "None" values (None = no BIN assignment)
        $self->bins( [ keys %{{ map { $_ => 1 } grep { defined $_ && $_ ne '' && $_ ne 'None' && $_ =~ /^BOLD:/ } map { $_->bin_uri } @{ $self->records } }} ] );
        
        # Count records with vs without valid BINs for logging
        my $total_records = $self->n_records;
        my $records_with_bins = scalar(grep { defined $_->bin_uri && $_->bin_uri ne '' && $_->bin_uri ne 'None' && $_->bin_uri =~ /^BOLD:/ } @{ $self->records });
        my $records_without_bins = $total_records - $records_with_bins;
        
        $log->info("Found " . $self->n_bins . " distinct BINs for $name");
        if ($records_without_bins > 0) {
            $log->info("Note: $records_without_bins of $total_records records lack valid BIN assignments");
        }
        $log->info("Found " . $self->n_bins . " distinct BINs for $name");
    }
    return $self->{'taxon'};
}

sub n_records { scalar @{ shift->records } }

sub n_bins { scalar @{ shift->bins } }

sub grade {
    my $self = shift;

    # If this list's size is greater than zero, then the taxon shares a BIN with another species
    my $is_shared = scalar($self->sharing_taxa);

    # Grade A means: >=10 specimens, in 1 unshared BIN
    if ( $self->n_records >= 10 && $self->n_bins == 1 && !$is_shared ) {
        $log->info($self->taxon->name . " is BAGS grade A");
        return 'A';
    }

    # Grade B means: 3-9 specimens, in 1 unshared BIN
    elsif ( 3 <= $self->n_records && $self->n_records < 10 && $self->n_bins == 1 && !$is_shared ) {
        $log->info($self->taxon->name . " is BAGS grade B");
        return 'B';
    }

    # Grade C means: more than 1 unshared BIN
    elsif ( $self->n_bins > 1 && !$is_shared ) {
        $log->info($self->taxon->name . " is BAGS grade C");
        return 'C';
    }

    # Grade D means: <3 specimens, in 1 unshared BIN
    elsif ( $self->n_records < 3 && $self->n_bins == 1 && !$is_shared ) {
        $log->info($self->taxon->name . " is BAGS grade D");
        return 'D';
    }

    # Grade E means: BIN sharing with other species
    elsif ($is_shared) {
        $log->info($self->taxon->name . " is BAGS grade E");
        return 'E';
    }

    # Grade F means: no valid BINs (all records lack BIN assignments)
    elsif ( $self->n_bins == 0 ) {
        $log->info($self->taxon->name . " is BAGS grade F (no valid BINs)");
        return 'F';
    }

    # This shouldn't happen
    else {
        $log->warn("Could not determine BAGS grade for " . $self->taxon->name . 
                  " (records: " . $self->n_records . ", bins: " . $self->n_bins . ", shared: $is_shared)");
        return 'F';
    }
}

sub taxa_sharing_bin {
    my ( $self, $bin ) = @_;
    my $orm = $self->taxon->result_source->schema;

    # Get all records with the same BIN URI that are identified to species level
    # We need to join with taxa table to ensure we only get species-level identifications
    my $bin_records = $orm->resultset('Bold')->search({ 
        bin_uri => $bin 
    });
    my $total_count = $bin_records->count;
    $log->info("Assessing $total_count total records sharing bin $bin");

    # Get all distinct species-level taxa for the records matching the BIN URI
    my %seen_species;
    my $processed_count = 0;
    my $valid_taxonid_count = 0;
    my $species_level_count = 0;
    my $higher_level_count = 0;
    
    while ( my $record = $bin_records->next ) {
        $processed_count++;
        
        # Get the taxon for this record - with proper error handling
        my $record_taxonid = $record->get_column('taxonid');  # Use get_column to get raw value
        
        # DEBUG: Enhanced logging for troubleshooting
        if ($processed_count <= 10) {
            $log->info("DEBUG Record $processed_count:");
            $log->info("  recordid: " . ($record->recordid // 'undef'));
            $log->info("  taxonid: " . ($record_taxonid // 'undef')); 
            $log->info("  bin_uri: " . ($record->bin_uri // 'undef'));
        }
        
        # FIXED: More robust taxonid validation
        # Handle both string and numeric taxonids, and strip whitespace
        my $clean_taxonid;
        if (defined $record_taxonid) {
            $clean_taxonid = "$record_taxonid";  # Stringify
            $clean_taxonid =~ s/^\s+|\s+$//g;    # Strip whitespace
        }
        
        unless (defined $clean_taxonid && $clean_taxonid ne '' && $clean_taxonid =~ /^\d+$/ && $clean_taxonid > 0) {
            if ($processed_count <= 10) {
                $log->info("  -> SKIP: Invalid taxonid (original: " . ($record_taxonid // 'undef') . 
                          ", cleaned: " . ($clean_taxonid // 'undef') . ")");
            }
            next;
        }
        $valid_taxonid_count++;
        
        if ($processed_count <= 10) {
            $log->info("  -> Valid taxonid: $clean_taxonid");
        }
        
        my $record_taxon;
        eval {
            # Use the cleaned taxonid for lookup
            my $taxon_rs = $orm->resultset('Taxa')->search({ taxonid => $clean_taxonid });
            $record_taxon = $taxon_rs->first;
        };
        if ($@) {
            $log->warn("  -> ERROR: Cannot find taxon for taxonid $clean_taxonid: $@");
            next;
        }
        
        if ($record_taxon) {
            if ($processed_count <= 10) {
                $log->info("  -> Found taxon: " . $record_taxon->name . " (level: " . $record_taxon->level . ")");
            }
            
            # Count records at species level AND subspecies level that could represent 
            # additional species diversity
            if ( $record_taxon->level eq 'species' ) {
                $species_level_count++;
                my $species_name = $record_taxon->name;
                $seen_species{$species_name}++;
                if ($processed_count <= 10) {
                    $log->info("  -> SPECIES: $species_name");
                }
            } elsif ( $record_taxon->level eq 'subspecies' ) {
                # For subspecies, extract the species name (first two words)
                # e.g., "Argia fumipennis violacea" -> "Argia fumipennis"
                my $subspecies_full_name = $record_taxon->name;
                my @name_parts = split /\s+/, $subspecies_full_name;
                if (@name_parts >= 2) {
                    my $parent_species_name = join(' ', @name_parts[0,1]);
                    $species_level_count++;  # Count subspecies as species-level diversity
                    $seen_species{$parent_species_name}++;
                    if ($processed_count <= 10) {
                        $log->info("  -> SUBSPECIES: $subspecies_full_name -> counting as $parent_species_name");
                    }
                } else {
                    if ($processed_count <= 10) {
                        $log->info("  -> SUBSPECIES: $subspecies_full_name (cannot parse parent species)");
                    }
                }
            } else {
                # For higher-level assignments, we need to check if this represents
                # additional species diversity within the BIN
                $higher_level_count++;
                if ($processed_count <= 10) {
                    $log->info("  -> HIGHER LEVEL: " . $record_taxon->level . " = " . $record_taxon->name);
                }
                # Note: We don't count higher-level assignments as separate species
                # because they don't represent confirmed species identifications
            }
        } else {
            if ($processed_count <= 10) {
                $log->warn("  -> ERROR: No taxon found for taxonid $clean_taxonid");
            }
        }
    }
    
    $log->info("SUMMARY: Processed $processed_count records");
    $log->info("  Valid taxonids: $valid_taxonid_count"); 
    $log->info("  Species + Subspecies level: $species_level_count");
    $log->info("  Higher level: $higher_level_count");
    $log->info("  Distinct species found: " . scalar(keys %seen_species));

    # Remove the current species
    my $current_species = $self->taxon->name;
    delete $seen_species{$current_species};

    # Report and return the names
    my @names = keys %seen_species;
    if (@names) {
        $log->info("Found " . scalar(@names) . " other species sharing bin $bin: " . 
                   join(", ", sort @names));
    } else {
        $log->info("Found 0 other species sharing bin $bin");
        # If we have higher-level assignments but no species-level ones,
        # this might indicate cryptic diversity that should be flagged
        if ($higher_level_count > 0 && $species_level_count == 0) {
            $log->warn("WARNING: BIN $bin has $higher_level_count higher-level assignments " .
                      "but no species-level identifications - possible cryptic diversity");
        }
    }
    
    return @names;
}

sub sharing_taxa {
    my $self = shift;
    my $orm  = $self->taxon->result_source->schema;

    # Get all species that share any BIN with this species
    my %all_sharing_species;
    
    for my $bin ( @{ $self->bins } ) {
        next unless defined $bin;
        my @sharing = $self->taxa_sharing_bin($bin);
        for my $species (@sharing) {
            $all_sharing_species{$species}++;
        }
    }

    # Report and return distinct species names
    my @names = keys %all_sharing_species;
    $log->info("Found " . scalar(@names) . " total species sharing " . $self->n_bins . " distinct BINs with " . $self->taxon->name);
    return @names;
}

sub AUTOLOAD {
    my ( $self, $arg ) = @_;
    my $method = $AUTOLOAD; # Contains the fully qualified name of the called method
    $method =~ s/.*:://;    # Remove package name

    # Check if the name exists in the object's hash as defined in the constructor
    if ( exists $self->{$method} ) {

        # If an argument is passed, set the value of the key in the object's hash
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
