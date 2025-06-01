#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use DBI;
use Time::HiRes qw(time);

# Haplotype identification within BINs and species groups
# Handles large datasets efficiently by processing in chunks

my $db_file;
my $log_level = 'INFO';
my $min_overlap = 0.6;  # 60% minimum overlap
my $chunk_size = 1000;  # Process in chunks for memory efficiency

GetOptions(
    'db=s'          => \$db_file,
    'log=s'         => \$log_level,
    'min-overlap=f' => \$min_overlap,
    'chunk-size=i'  => \$chunk_size,
);

die "Database file required (--db)\n" unless $db_file;

sub log_msg {
    my ($level, $msg) = @_;
    return if $log_level eq 'ERROR' && $level ne 'ERROR';
    return if $log_level eq 'WARN' && $level eq 'INFO';
    
    my $timestamp = scalar(localtime);
    print STDERR "[$timestamp] $level: $msg\n";
}

# Connect to database
log_msg('INFO', "Connecting to database: $db_file");
my $dbh = DBI->connect("dbi:SQLite:$db_file", "", "", {
    RaiseError => 1,
    AutoCommit => 1,  # Temporarily enable AutoCommit for PRAGMA statements
    sqlite_see_if_its_a_number => 1,
}) or die "Cannot connect to database: $DBI::errstr";

# Enable performance optimizations (must be done outside of transactions)
$dbh->do("PRAGMA journal_mode = WAL");
$dbh->do("PRAGMA synchronous = NORMAL");
$dbh->do("PRAGMA cache_size = -256000");  # 256MB cache
$dbh->do("PRAGMA temp_store = MEMORY");

# Now disable AutoCommit for transaction control
$dbh->{AutoCommit} = 0;

# Create haplotype tables if they don't exist
create_haplotype_tables($dbh);

# Get sequences grouped by BIN and species
log_msg('INFO', "Retrieving sequences for haplotype analysis...");

# First process BIN-based groups
my $bin_sql = q{
    SELECT recordid, bin_uri, species, nuc, nuc_basecount
    FROM bold 
    WHERE bin_uri IS NOT NULL 
      AND bin_uri != 'None' 
      AND bin_uri != ''
      AND nuc IS NOT NULL 
      AND nuc != ''
      AND species IS NOT NULL
      AND species != ''
    ORDER BY bin_uri, species
};

# Then process species-based groups (no BIN)
my $species_sql = q{
    SELECT recordid, species, nuc, nuc_basecount
    FROM bold 
    WHERE (bin_uri IS NULL OR bin_uri = 'None' OR bin_uri = '')
      AND nuc IS NOT NULL 
      AND nuc != ''
      AND species IS NOT NULL
      AND species != ''
    ORDER BY species
};

# Process BIN-based haplotypes first
log_msg('INFO', "Processing BIN-based haplotypes...");
my %global_haplotypes = ();  # Track all haplotypes globally for deduplication
process_bin_haplotypes($dbh, \%global_haplotypes);

# Then process species-based haplotypes
log_msg('INFO', "Processing species-based haplotypes...");
process_species_haplotypes($dbh, \%global_haplotypes);

# Generate assessment output
log_msg('INFO', "Generating haplotype assessment output...");
generate_assessment_output($dbh);

$dbh->commit;
$dbh->disconnect;

log_msg('INFO', "Haplotype analysis completed successfully");

# Subroutines

sub create_haplotype_tables {
    my ($dbh) = @_;
    
    log_msg('INFO', "Creating haplotype tables...");
    
    $dbh->do(q{
        CREATE TABLE IF NOT EXISTS "haplotypes" (
            "haplotype_id" INTEGER PRIMARY KEY,
            "haplotype_name" TEXT NOT NULL UNIQUE,
            "bin_uri" TEXT,
            "species_name" TEXT,
            "sequence_length" INTEGER,
            "record_count" INTEGER DEFAULT 1,
            "created_date" TIMESTAMP DEFAULT CURRENT_TIMESTAMP
        )
    });
    
    $dbh->do(q{
        CREATE TABLE IF NOT EXISTS "bold_haplotypes" (
            "bold_haplotype_id" INTEGER PRIMARY KEY,
            "recordid" INTEGER NOT NULL,
            "haplotype_id" INTEGER NOT NULL,
            "is_representative" INTEGER DEFAULT 0,
            FOREIGN KEY(recordid) REFERENCES bold(recordid),
            FOREIGN KEY(haplotype_id) REFERENCES haplotypes(haplotype_id),
            CONSTRAINT unique_record_haplotype UNIQUE (recordid)
        )
    });
    
    # Clear existing haplotype data for fresh analysis
    $dbh->do("DELETE FROM bold_haplotypes");
    $dbh->do("DELETE FROM haplotypes");
}

sub process_bin_haplotypes {
    my ($dbh, $global_haplotypes) = @_;
    
    my $sth = $dbh->prepare(q{
        SELECT recordid, bin_uri, species, nuc, nuc_basecount
        FROM bold 
        WHERE bin_uri IS NOT NULL 
          AND bin_uri != 'None' 
          AND bin_uri != ''
          AND nuc IS NOT NULL 
          AND nuc != ''
          AND species IS NOT NULL
          AND species != ''
        ORDER BY bin_uri, species
    });
    
    $sth->execute();
    
    my %current_group = ();
    my $current_bin = '';
    my $current_species = '';
    my $processed_records = 0;
    
    while (my $row = $sth->fetchrow_hashref()) {
        my $bin_species_key = "$row->{bin_uri}|$row->{species}";
        
        # If we've moved to a new BIN+species group, process the previous group
        if ($current_bin ne $row->{bin_uri} || $current_species ne $row->{species}) {
            if (%current_group) {
                process_haplotype_group(\%current_group, $current_bin, $current_species, $dbh, $global_haplotypes);
                %current_group = ();
            }
            $current_bin = $row->{bin_uri};
            $current_species = $row->{species};
        }
        
        # Add record to current group
        $current_group{$row->{recordid}} = {
            sequence => normalize_sequence($row->{nuc}),
            original_seq => $row->{nuc},
            length => $row->{nuc_basecount} || length($row->{nuc}),
        };
        
        $processed_records++;
        if ($processed_records % 10000 == 0) {
            log_msg('INFO', "Processed $processed_records BIN-based records...");
        }
    }
    
    # Process final group
    if (%current_group) {
        process_haplotype_group(\%current_group, $current_bin, $current_species, $dbh, $global_haplotypes);
    }
    
    log_msg('INFO', "Completed BIN-based processing: $processed_records records");
}

sub process_species_haplotypes {
    my ($dbh, $global_haplotypes) = @_;
    
    my $sth = $dbh->prepare(q{
        SELECT recordid, species, nuc, nuc_basecount
        FROM bold 
        WHERE (bin_uri IS NULL OR bin_uri = 'None' OR bin_uri = '')
          AND nuc IS NOT NULL 
          AND nuc != ''
          AND species IS NOT NULL
          AND species != ''
        ORDER BY species
    });
    
    $sth->execute();
    
    my %current_group = ();
    my $current_species = '';
    my $processed_records = 0;
    
    while (my $row = $sth->fetchrow_hashref()) {
        # If we've moved to a new species, process the previous group
        if ($current_species ne $row->{species}) {
            if (%current_group) {
                process_haplotype_group(\%current_group, undef, $current_species, $dbh, $global_haplotypes);
                %current_group = ();
            }
            $current_species = $row->{species};
        }
        
        # Add record to current group
        $current_group{$row->{recordid}} = {
            sequence => normalize_sequence($row->{nuc}),
            original_seq => $row->{nuc},
            length => $row->{nuc_basecount} || length($row->{nuc}),
        };
        
        $processed_records++;
        if ($processed_records % 10000 == 0) {
            log_msg('INFO', "Processed $processed_records species-based records...");
        }
    }
    
    # Process final group
    if (%current_group) {
        process_haplotype_group(\%current_group, undef, $current_species, $dbh, $global_haplotypes);
    }
    
    log_msg('INFO', "Completed species-based processing: $processed_records records");
}

sub process_haplotype_group {
    my ($group, $bin_uri, $species, $dbh, $global_haplotypes) = @_;
    
    my @record_ids = keys %$group;
    return unless @record_ids;
    
    log_msg('INFO', sprintf("Processing group: %s (%d records)", 
        $bin_uri ? "$bin_uri|$species" : $species, 
        scalar @record_ids
    ));
    
    my @haplotype_groups = ();
    
    # Group sequences into haplotypes
    for my $record_id (@record_ids) {
        my $record = $group->{$record_id};
        my $assigned = 0;
        
        # Try to assign to existing haplotype group
        for my $hap_group (@haplotype_groups) {
            if (sequences_match($record->{sequence}, $hap_group->[0]{sequence})) {
                push @$hap_group, {
                    record_id => $record_id,
                    sequence => $record->{sequence},
                    original_seq => $record->{original_seq},
                    length => $record->{length},
                };
                $assigned = 1;
                last;
            }
        }
        
        # Create new haplotype group if not assigned
        unless ($assigned) {
            push @haplotype_groups, [{
                record_id => $record_id,
                sequence => $record->{sequence},
                original_seq => $record->{original_seq},
                length => $record->{length},
            }];
        }
    }
    
    # Create haplotype records and assignments
    my $haplotype_counter = 1;
    for my $hap_group (@haplotype_groups) {
        my $haplotype_name = generate_haplotype_name($bin_uri, $species, $haplotype_counter);
        
        # Check for global deduplication (for species-based haplotypes)
        my $existing_haplotype_id = check_global_haplotype($hap_group->[0]{sequence}, $species, $global_haplotypes);
        
        my $haplotype_id;
        if ($existing_haplotype_id) {
            $haplotype_id = $existing_haplotype_id;
            log_msg('INFO', "Reusing existing haplotype for $haplotype_name");
        } else {
            $haplotype_id = create_haplotype($dbh, $haplotype_name, $bin_uri, $species, $hap_group);
            
            # Add to global tracking
            $global_haplotypes->{$hap_group->[0]{sequence}} = {
                haplotype_id => $haplotype_id,
                species => $species,
                bin_uri => $bin_uri,
            };
        }
        
        # Assign all records in this group to the haplotype
        assign_records_to_haplotype($dbh, $hap_group, $haplotype_id);
        
        $haplotype_counter++;
    }
}

sub normalize_sequence {
    my ($seq) = @_;
    return '' unless defined $seq;
    
    # Remove whitespace and convert to uppercase
    $seq =~ s/\s+//g;
    $seq = uc($seq);
    
    # Keep only valid IUPAC nucleotide codes (including ambiguous bases)
    # A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N
    $seq =~ s/[^ACGTRYSWKMBDHVNU]//g;
    
    # Convert U to T for consistency
    $seq =~ tr/U/T/;
    
    return $seq;
}

sub reverse_complement {
    my ($seq) = @_;
    $seq =~ tr/ATCG/TAGC/;
    return reverse($seq);
}

sub sequences_match {
    my ($seq1, $seq2) = @_;
    
    return 0 unless $seq1 && $seq2;
    
    # Try both orientations
    return 1 if sequences_overlap($seq1, $seq2);
    return 1 if sequences_overlap($seq1, reverse_complement($seq2));
    
    return 0;
}

sub sequences_overlap {
    my ($seq1, $seq2) = @_;
    
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    
    return 0 if $len1 == 0 || $len2 == 0;
    
    my $min_len = $len1 < $len2 ? $len1 : $len2;
    my $required_overlap = int($min_len * $min_overlap);
    
    return 0 if $required_overlap < 1;
    
    # Check if sequences are identical in their overlap region
    if ($len1 <= $len2) {
        # seq1 is shorter or equal, check if it matches the beginning of seq2
        return substr($seq2, 0, $len1) eq $seq1;
    } else {
        # seq2 is shorter, check if it matches the beginning of seq1
        return substr($seq1, 0, $len2) eq $seq2;
    }
}

sub check_global_haplotype {
    my ($sequence, $species, $global_haplotypes) = @_;
    
    # Only check for species-based sequences (no BIN)
    for my $existing_seq (keys %$global_haplotypes) {
        my $existing = $global_haplotypes->{$existing_seq};
        
        # Only match against BIN-based haplotypes of the same species
        if ($existing->{species} eq $species && 
            $existing->{bin_uri} && 
            sequences_match($sequence, $existing_seq)) {
            return $existing->{haplotype_id};
        }
    }
    
    return undef;
}

sub generate_haplotype_name {
    my ($bin_uri, $species, $counter) = @_;
    
    if ($bin_uri) {
        return "${bin_uri}_H${counter}";
    } else {
        # Clean species name for use in identifier
        my $clean_species = $species;
        $clean_species =~ s/\s+/_/g;
        $clean_species =~ s/[^A-Za-z0-9_]//g;
        return "${clean_species}_H${counter}";
    }
}

sub create_haplotype {
    my ($dbh, $name, $bin_uri, $species, $hap_group) = @_;
    
    my $representative = $hap_group->[0];  # Use first record as representative
    
    my $sth = $dbh->prepare(q{
        INSERT INTO haplotypes (haplotype_name, bin_uri, species_name, sequence_length, record_count)
        VALUES (?, ?, ?, ?, ?)
    });
    
    $sth->execute($name, $bin_uri, $species, $representative->{length}, scalar @$hap_group);
    
    return $dbh->last_insert_id("", "", "haplotypes", "haplotype_id");
}

sub assign_records_to_haplotype {
    my ($dbh, $hap_group, $haplotype_id) = @_;
    
    my $sth = $dbh->prepare(q{
        INSERT INTO bold_haplotypes (recordid, haplotype_id, is_representative)
        VALUES (?, ?, ?)
    });
    
    for my $i (0 .. $#{$hap_group}) {
        my $record = $hap_group->[$i];
        my $is_representative = ($i == 0) ? 1 : 0;  # First record is representative
        
        $sth->execute($record->{record_id}, $haplotype_id, $is_representative);
    }
}

sub generate_assessment_output {
    my ($dbh) = @_;
    
    # Output header compatible with other assessment scripts
    print "recordid\thaplotype_id\tstatus\tnotes\n";
    
    my $sth = $dbh->prepare(q{
        SELECT bh.recordid, h.haplotype_name
        FROM bold_haplotypes bh
        JOIN haplotypes h ON bh.haplotype_id = h.haplotype_id
        ORDER BY bh.recordid
    });
    
    $sth->execute();
    
    while (my ($recordid, $haplotype_name) = $sth->fetchrow_array()) {
        print "$recordid\t$haplotype_name\t1\tHaplotype assigned\n";
    }
}