#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use DBI;
use Time::HiRes qw(time);
use Digest::MD5 qw(md5_hex);
use Parallel::ForkManager;

# Optimized Haplotype identification within BINs and species groups
# Major optimizations:
# 1. Hash-based sequence deduplication for exact matches
# 2. Parallel processing of taxonomic groups
# 3. Chunked database processing with better indexing
# 4. Optimized sequence comparison algorithms
# 5. Memory-efficient processing

my $db_file;
my $log_level = 'INFO';
my $min_overlap = 0.6;  # 60% minimum overlap
my $chunk_size = 5000;  # Increased chunk size for efficiency
my $max_processes = 4;  # Parallel processing (adjust based on your system)
my $hash_cache_size = 100000;  # Hash cache for frequent sequences

GetOptions(
    'db=s'           => \$db_file,
    'log=s'          => \$log_level,
    'min-overlap=f'  => \$min_overlap,
    'chunk-size=i'   => \$chunk_size,
    'processes=i'    => \$max_processes,
    'cache-size=i'   => \$hash_cache_size,
);

die "Database file required (--db)\n" unless $db_file;

# Global data structures with optimization
my %sequence_hashes = ();  # sequence_hash -> canonical_sequence
my %hash_to_haplotype = ();  # sequence_hash -> haplotype_name
my %haplotype_assignments = ();  # recordid -> haplotype_name
my %haplotype_names_used = ();  # Track used names
my %reverse_complement_cache = ();  # Cache for reverse complements

sub log_msg {
    my ($level, $msg) = @_;
    return if $log_level eq 'ERROR' && $level ne 'ERROR';
    return if $log_level eq 'WARN' && $level eq 'INFO';
    
    my $timestamp = scalar(localtime);
    print STDERR "[$timestamp] $level: $msg\n";
}

# Connect to database with optimizations
log_msg('INFO', "Connecting to database: $db_file");
my $dbh = DBI->connect("dbi:SQLite:$db_file", "", "", {
    RaiseError => 1,
    AutoCommit => 1,
    ReadOnly => 1,
    sqlite_see_if_its_a_number => 1,
    sqlite_cache_size => 100000,  # Increase cache size
    sqlite_temp_store => 2,       # Use memory for temp storage
}) or die "Cannot connect to database: $DBI::errstr";

# Apply query optimizations
$dbh->do("PRAGMA query_planner = 1");
$dbh->do("PRAGMA optimize");

log_msg('INFO', "Starting optimized haplotype analysis...");
my $start_time = time();

# Get total record count for progress tracking
my $total_count = $dbh->selectrow_array(q{
    SELECT COUNT(*) FROM bold 
    WHERE nuc IS NOT NULL AND nuc != '' 
    AND species IS NOT NULL AND species != ''
});

log_msg('INFO', "Total records to process: $total_count");

# Process BIN-based groups with optimization
log_msg('INFO', "Processing BIN-based haplotypes with parallel processing...");
process_bin_haplotypes_optimized($dbh);

# Process species-based groups  
log_msg('INFO', "Processing species-based haplotypes...");
process_species_haplotypes_optimized($dbh);

# Generate TSV output
log_msg('INFO', "Generating haplotype assessment output...");
generate_tsv_output(\%haplotype_assignments);

$dbh->disconnect;

my $end_time = time();
my $duration = $end_time - $start_time;
log_msg('INFO', sprintf("Haplotype analysis completed in %.2f seconds", $duration));
log_msg('INFO', sprintf("Processing rate: %.0f records/second", $total_count / $duration));

# Optimized subroutines

sub process_bin_haplotypes_optimized {
    my ($dbh) = @_;
    
    # Get distinct BIN+species combinations for parallel processing
    my $sth = $dbh->prepare(q{
        SELECT bin_uri, species, COUNT(*) as record_count
        FROM bold 
        WHERE bin_uri IS NOT NULL 
          AND bin_uri != 'None' 
          AND bin_uri != ''
          AND nuc IS NOT NULL 
          AND nuc != ''
          AND species IS NOT NULL
          AND species != ''
        GROUP BY bin_uri, species
        ORDER BY record_count DESC
    });
    
    $sth->execute();
    my @groups = ();
    
    while (my $row = $sth->fetchrow_hashref()) {
        push @groups, {
            bin_uri => $row->{bin_uri},
            species => $row->{species},
            count => $row->{record_count}
        };
    }
    
    log_msg('INFO', sprintf("Found %d BIN+species groups for processing", scalar @groups));
    
    # Process groups in parallel
    my $pm = Parallel::ForkManager->new($max_processes);
    
    # Data sharing between processes
    my %child_results = ();
    $pm->run_on_finish(sub {
        my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
        if (defined $data_ref) {
            # Merge results from child process
            my ($assignments, $names_used) = @$data_ref;
            %haplotype_assignments = (%haplotype_assignments, %$assignments);
            %haplotype_names_used = (%haplotype_names_used, %$names_used);
        }
    });
    
    my $processed_groups = 0;
    foreach my $group (@groups) {
        $pm->start and next;  # Fork child process
        
        # Child process
        my %local_assignments = ();
        my %local_names_used = ();
        
        process_single_group_optimized($dbh, $group, \%local_assignments, \%local_names_used);
        
        # Return results to parent
        $pm->finish(0, [\%local_assignments, \%local_names_used]);
    }
    
    $pm->wait_all_children;
    
    log_msg('INFO', sprintf("Completed BIN-based processing: %d groups", scalar @groups));
}

sub process_species_haplotypes_optimized {
    my ($dbh) = @_;
    
    # Get species without BINs for processing
    my $sth = $dbh->prepare(q{
        SELECT species, COUNT(*) as record_count
        FROM bold 
        WHERE (bin_uri IS NULL OR bin_uri = 'None' OR bin_uri = '')
          AND nuc IS NOT NULL 
          AND nuc != ''
          AND species IS NOT NULL
          AND species != ''
        GROUP BY species
        ORDER BY record_count DESC
    });
    
    $sth->execute();
    my @species_groups = ();
    
    while (my $row = $sth->fetchrow_hashref()) {
        push @species_groups, {
            species => $row->{species},
            count => $row->{record_count}
        };
    }
    
    log_msg('INFO', sprintf("Found %d species groups for processing", scalar @species_groups));
    
    # Process smaller groups sequentially for efficiency
    my $processed = 0;
    foreach my $group (@species_groups) {
        my %local_assignments = ();
        my %local_names_used = ();
        
        process_single_group_optimized($dbh, $group, \%local_assignments, \%local_names_used);
        
        # Merge results
        %haplotype_assignments = (%haplotype_assignments, %local_assignments);
        %haplotype_names_used = (%haplotype_names_used, %local_names_used);
        
        $processed++;
        if ($processed % 100 == 0) {
            log_msg('INFO', "Processed $processed species groups...");
        }
    }
    
    log_msg('INFO', sprintf("Completed species-based processing: %d groups", scalar @species_groups));
}

sub process_single_group_optimized {
    my ($dbh, $group, $assignments, $names_used) = @_;
    
    my $bin_uri = $group->{bin_uri} || undef;
    my $species = $group->{species};
    
    # Fetch records for this group with optimized query
    my $sql = $bin_uri ? 
        q{SELECT recordid, nuc, nuc_basecount FROM bold 
          WHERE bin_uri = ? AND species = ? 
          AND nuc IS NOT NULL AND nuc != ''} :
        q{SELECT recordid, nuc, nuc_basecount FROM bold 
          WHERE (bin_uri IS NULL OR bin_uri = 'None' OR bin_uri = '') 
          AND species = ? AND nuc IS NOT NULL AND nuc != ''};
    
    my $sth = $dbh->prepare($sql);
    if ($bin_uri) {
        $sth->execute($bin_uri, $species);
    } else {
        $sth->execute($species);
    }
    
    # Hash-based sequence clustering
    my %sequence_groups = ();  # hash -> [record_info, ...]
    my %processed_hashes = (); # Track already processed sequence hashes
    
    while (my $row = $sth->fetchrow_hashref()) {
        my $normalized_seq = normalize_sequence($row->{nuc});
        next unless $normalized_seq;
        
        my $seq_hash = md5_hex($normalized_seq);
        
        # Check if we've already processed this exact sequence
        if (exists $processed_hashes{$seq_hash}) {
            # Assign to existing haplotype
            $assignments->{$row->{recordid}} = $processed_hashes{$seq_hash};
            next;
        }
        
        # Check reverse complement hash
        my $rev_comp = get_reverse_complement_cached($normalized_seq);
        my $rev_hash = md5_hex($rev_comp);
        
        if (exists $processed_hashes{$rev_hash}) {
            # Assign to existing haplotype (reverse complement match)
            $assignments->{$row->{recordid}} = $processed_hashes{$rev_hash};
            $processed_hashes{$seq_hash} = $processed_hashes{$rev_hash};
            next;
        }
        
        # New sequence - check for similar sequences
        my $assigned_haplotype = find_matching_haplotype(
            $normalized_seq, \%sequence_groups, $assignments
        );
        
        if ($assigned_haplotype) {
            $assignments->{$row->{recordid}} = $assigned_haplotype;
            $processed_hashes{$seq_hash} = $assigned_haplotype;
        } else {
            # Create new haplotype
            my $new_haplotype = generate_unique_haplotype_name(
                $bin_uri, $species, scalar(keys %sequence_groups) + 1, $names_used
            );
            
            $sequence_groups{$seq_hash} = {
                haplotype_name => $new_haplotype,
                sequence => $normalized_seq,
                records => [$row->{recordid}]
            };
            
            $assignments->{$row->{recordid}} = $new_haplotype;
            $processed_hashes{$seq_hash} = $new_haplotype;
            $names_used->{$new_haplotype} = 1;
        }
    }
}

sub find_matching_haplotype {
    my ($sequence, $sequence_groups, $assignments) = @_;
    
    # Only check against representative sequences from each group
    for my $hash (keys %$sequence_groups) {
        my $group = $sequence_groups->{$hash};
        if (sequences_match_optimized($sequence, $group->{sequence})) {
            return $group->{haplotype_name};
        }
    }
    
    return undef;
}

sub normalize_sequence {
    my ($seq) = @_;
    return '' unless defined $seq;
    
    # Optimized normalization
    $seq =~ tr/a-z/A-Z/;           # Convert to uppercase
    $seq =~ s/\s+//g;              # Remove whitespace
    $seq =~ s/[^ACGTRYSWKMBDHVNU]//g;  # Keep only valid IUPAC codes
    $seq =~ tr/U/T/;               # Convert U to T
    
    return $seq;
}

sub get_reverse_complement_cached {
    my ($seq) = @_;
    
    return $reverse_complement_cache{$seq} if exists $reverse_complement_cache{$seq};
    
    my $rev_comp = reverse_complement($seq);
    
    # Limit cache size to prevent memory bloat
    if (keys %reverse_complement_cache > $hash_cache_size) {
        %reverse_complement_cache = ();
    }
    
    $reverse_complement_cache{$seq} = $rev_comp;
    return $rev_comp;
}

sub reverse_complement {
    my ($seq) = @_;
    $seq =~ tr/ATCG/TAGC/;
    return reverse($seq);
}

sub sequences_match_optimized {
    my ($seq1, $seq2) = @_;
    
    return 0 unless $seq1 && $seq2;
    
    # Quick length check
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    
    return 0 if $len1 == 0 || $len2 == 0;
    
    # If sequences are identical, return immediately
    return 1 if $seq1 eq $seq2;
    
    # Check reverse complement
    my $rev_comp = get_reverse_complement_cached($seq2);
    return 1 if $seq1 eq $rev_comp;
    
    # Check overlap for different length sequences
    return sequences_overlap_optimized($seq1, $seq2) || 
           sequences_overlap_optimized($seq1, $rev_comp);
}

sub sequences_overlap_optimized {
    my ($seq1, $seq2) = @_;
    
    my $len1 = length($seq1);
    my $len2 = length($seq2);
    
    return 0 if $len1 == 0 || $len2 == 0;
    
    my $min_len = $len1 < $len2 ? $len1 : $len2;
    my $required_overlap = int($min_len * $min_overlap);
    
    return 0 if $required_overlap < 1;
    
    # Optimized overlap check using substr
    if ($len1 <= $len2) {
        return substr($seq2, 0, $len1) eq $seq1;
    } else {
        return substr($seq1, 0, $len2) eq $seq2;
    }
}

sub generate_unique_haplotype_name {
    my ($bin_uri, $species, $counter, $names_used) = @_;
    
    my $base_name;
    if ($bin_uri) {
        $base_name = "${bin_uri}_H";
    } else {
        my $clean_species = $species;
        $clean_species =~ s/\s+/_/g;
        $clean_species =~ s/[^A-Za-z0-9_]//g;
        $base_name = "${clean_species}_H";
    }
    
    my $attempt_counter = $counter;
    my $candidate_name = "${base_name}${attempt_counter}";
    
    while (exists $names_used->{$candidate_name}) {
        $attempt_counter++;
        $candidate_name = "${base_name}${attempt_counter}";
        
        if ($attempt_counter > 10000) {
            log_msg('ERROR', "Unable to generate unique haplotype name for $base_name");
            last;
        }
    }
    
    return $candidate_name;
}

sub generate_tsv_output {
    my ($assignments) = @_;
    
    print "recordid\thaplotype_id\tstatus\tnotes\n";
    
    my $count = 0;
    for my $recordid (sort { $a <=> $b } keys %$assignments) {
        my $haplotype_name = $assignments->{$recordid};
        print "$recordid\t$haplotype_name\t1\tHaplotype assigned (optimized)\n";
        $count++;
    }
    
    log_msg('INFO', "Generated $count haplotype assignments");
}
