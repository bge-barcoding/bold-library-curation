use strict;
use warnings;
use BCDM::ORM;
use Getopt::Long;
use Log::Log4perl qw(:easy);
use DBI;

# Levels in the BOLD taxonomy
my @levels = qw[
    kingdom
    phylum
    class
    order
    family
    subfamily
    genus
    species
    subspecies
];

# Process command line arguments
my $db_file; # where to access database file
my $log_level = 'INFO'; # verbosity level for logger
my $chunk_size = 10000; # records to process in each chunk
GetOptions(
    'db=s'  => \$db_file,
    'log=s' => \$log_level,
    'chunk=i' => \$chunk_size,
);

# Initialize Log::Log4perl
Log::Log4perl->init(\<<"END");
  log4perl.rootLogger = $log_level, Screen
  log4perl.appender.Screen = Log::Log4perl::Appender::Screen
  log4perl.appender.Screen.stderr = 1
  log4perl.appender.Screen.layout = Log::Log4perl::Layout::SimpleLayout
END

# Instantiate logger
my $log = Log::Log4perl->get_logger('load_taxonomy');

# Connect to the database with optimized settings
$log->info("Going to connect to database $db_file");
my $schema = BCDM::ORM->connect("dbi:SQLite:$db_file", "", "", { 
    quote_char => '"',
    sqlite_use_immediate_transaction => 1,
});

# Get direct DBI handle for raw SQL operations
my $dbh = $schema->storage->dbh;

# Optimize SQLite for bulk operations
$log->info("Optimizing database for bulk operations");
$dbh->do("PRAGMA journal_mode = WAL");
$dbh->do("PRAGMA synchronous = OFF");
$dbh->do("PRAGMA cache_size = 100000");
$dbh->do("PRAGMA temp_store = MEMORY");

# PASS 1: Extract all unique taxonomic paths and bulk insert them
$log->info("PASS 1: Extracting and creating unique taxonomic paths");

my %unique_paths; # Move to global scope
extract_unique_taxonomic_paths(\%unique_paths);
$log->info("Found " . scalar(keys %unique_paths) . " unique taxonomic paths");

bulk_insert_taxonomy(\%unique_paths);

# PASS 2: Bulk update BOLD records with taxonomy IDs
$log->info("PASS 2: Linking BOLD records to taxonomy");
bulk_update_bold_records();

# Restore normal SQLite settings
$log->info("Restoring normal database settings");
$dbh->do("PRAGMA synchronous = FULL");

$log->info("Taxonomy loading completed");

sub extract_unique_taxonomic_paths {
    my $unique_paths = shift; # Accept hash reference as parameter
    
    $log->info("Extracting unique taxonomic paths from BOLD records");
    
    # DEBUG: Check what columns are available
    my $sample_rs = $schema->resultset('Bold')->search({}, { rows => 1 });
    if (my $sample_record = $sample_rs->next) {
        my @available_columns = $sample_record->result_source->columns;
        $log->info("DEBUG: Available columns: " . join(", ", @available_columns));
        
        # Check if our expected columns exist
        for my $col ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species') {
            my $has_col = grep { $_ eq $col } @available_columns;
            $log->info("DEBUG: Column '$col' " . ($has_col ? "EXISTS" : "MISSING"));
        }
    }
    
    # Get all BOLD records in chunks to avoid memory issues
    my $total_records = $schema->resultset('Bold')->count;
    $log->info("Processing $total_records total records");
    
    my $offset = 0;
    while ($offset < $total_records) {
        my $rs = $schema->resultset('Bold')->search(
            {},
            {
                rows => $chunk_size,
                offset => $offset,
                columns => ['recordid', 'kingdom', @levels]
            }
        );
        
        my $record_count = 0;
        while (my $record = $rs->next) {
            $record_count++;
            my $kingdom = $record->get_column('kingdom') || '';
            
            # DEBUG: Log first few records
            if ($record_count <= 3) {
                $log->info("DEBUG Record $record_count: kingdom='$kingdom'");
                for my $level (@levels) {
                    my $value = $record->get_column($level) || '';
                    $log->info("DEBUG   $level='$value'");
                }
            }
            
            # Skip if no kingdom
            unless ($kingdom) {
                $log->info("DEBUG: Skipping record $record_count - no kingdom") if $record_count <= 10;
                next;
            }
            
            # Build taxonomic path from kingdom down to lowest defined level
            # FIXED: Collect all defined levels, skipping "None" values
            my $path = [];
            my $last_valid_rank = -1;
            
            for my $i (0 .. $#levels) {
                my $level = $levels[$i];
                my $name = $record->get_column($level) || '';
                
                # Skip empty, "None", or undefined values but continue processing
                if ($name eq '' || $name eq 'None' || !defined($name)) {
                    $log->info("DEBUG Record $record_count: Skipping empty/None $level") if $record_count <= 3;
                    next;
                }
                
                push @$path, {
                    kingdom => $kingdom,
                    name => $name,
                    level => $level,
                    rank => $i
                };
                $last_valid_rank = $i;
            }
            
            # DEBUG: Log path length for first few records
            if ($record_count <= 3) {
                $log->info("DEBUG Record $record_count: path length = " . scalar(@$path));
                for my $node (@$path) {
                    $log->info("DEBUG   Found: $node->{level} = $node->{name}");
                }
            }
            
            # Store each node in the path
            for my $i (0 .. $#$path) {
                my $node = $path->[$i];
                my $parent_key = $i > 0 ? create_path_key($path->[$i-1]) : 'ROOT';
                my $node_key = create_path_key($node);
                
                $unique_paths->{$node_key} = {
                    %$node,
                    parent_key => $parent_key
                };
            }
        }
        
        $offset += $chunk_size;
        $log->info("Processed $offset records, found " . scalar(keys %$unique_paths) . " unique paths so far") if $offset % ($chunk_size * 10) == 0;
    }
    
    $log->info("DEBUG: Final unique_paths count = " . scalar(keys %$unique_paths));
    if (scalar(keys %$unique_paths) > 0) {
        my @sample_keys = (keys %$unique_paths)[0..2];
        for my $key (@sample_keys) {
            $log->info("DEBUG: Sample path key: $key");
        }
    }
    
    # No return statement needed - we modified the hash reference directly
}

sub create_path_key {
    my $node = shift;
    return join('|', $node->{kingdom}, $node->{level}, $node->{name});
}

sub bulk_insert_taxonomy {
    my $unique_paths = shift;
    
    $log->info("Bulk inserting taxonomy nodes");
    
    # Sort by rank to ensure parents are inserted before children
    my @sorted_paths = sort { 
        $unique_paths->{$a}->{rank} <=> $unique_paths->{$b}->{rank} 
    } keys %$unique_paths;
    
    # Prepare statements
    my $find_sth = $dbh->prepare(
        "SELECT taxonid FROM taxa WHERE kingdom = ? AND name = ? AND level = ? AND " .
        "(parent_taxonid = ? OR (parent_taxonid IS NULL AND ? IS NULL))"
    );
    
    my $insert_sth = $dbh->prepare(
        "INSERT INTO taxa (kingdom, name, level, parent_taxonid) VALUES (?, ?, ?, ?)"
    );
    
    my %taxon_cache; # Maps path_key to taxonid
    my $created_count = 0;
    my $batch_count = 0;
    
    $dbh->begin_work;
    
    for my $path_key (@sorted_paths) {
        my $node = $unique_paths->{$path_key};
        my $parent_taxonid = undef;
        
        # Get parent taxonid if not root
        if ($node->{parent_key} ne 'ROOT') {
            $parent_taxonid = $taxon_cache{$node->{parent_key}};
        }
        
        # Check if taxon already exists
        $find_sth->execute(
            $node->{kingdom}, 
            $node->{name}, 
            $node->{level}, 
            $parent_taxonid,
            $parent_taxonid
        );
        
        my ($existing_taxonid) = $find_sth->fetchrow_array;
        
        if ($existing_taxonid) {
            $taxon_cache{$path_key} = $existing_taxonid;
        } else {
            # Insert new taxon
            $insert_sth->execute(
                $node->{kingdom},
                $node->{name},
                $node->{level},
                $parent_taxonid
            );
            
            my $new_taxonid = $dbh->last_insert_id("", "", "taxa", "taxonid");
            $taxon_cache{$path_key} = $new_taxonid;
            $created_count++;
        }
        
        $batch_count++;
        if ($batch_count % 1000 == 0) {
            $dbh->commit;
            $dbh->begin_work;
            $log->info("Processed $batch_count taxonomy nodes, created $created_count new ones");
        }
    }
    
    $dbh->commit;
    $log->info("Completed taxonomy insertion: $created_count new taxa created from " . 
               scalar(@sorted_paths) . " unique paths");
}

sub bulk_update_bold_records {
    $log->info("Bulk updating BOLD records with taxonomy IDs");
    
    # Prepare statements for finding taxonids and updating records
    my $update_sth = $dbh->prepare("UPDATE bold SET taxonid = ? WHERE recordid = ?");
    
    my $total_records = $schema->resultset('Bold')->count;
    my $offset = 0;
    my $updated_count = 0;
    
    while ($offset < $total_records) {
        my @updates;
        
        # Get chunk of BOLD records
        my $rs = $schema->resultset('Bold')->search(
            {},
            {
                rows => $chunk_size,
                offset => $offset,
                columns => ['recordid', 'kingdom', @levels]
            }
        );
        
        while (my $record = $rs->next) {
            my $kingdom = $record->get_column('kingdom') || '';
            next unless $kingdom;
            
            # Find the lowest defined taxon for this record
            my $taxonid = find_taxonid_for_record($record);
            
            if ($taxonid) {
                push @updates, [$taxonid, $record->get_column('recordid')];
            }
        }
        
        # Batch update
        if (@updates) {
            $dbh->begin_work;
            for my $update (@updates) {
                $update_sth->execute(@$update);
            }
            $dbh->commit;
            $updated_count += scalar(@updates);
        }
        
        $offset += $chunk_size;
        $log->info("Updated $updated_count of $total_records records") 
            if $offset % ($chunk_size * 10) == 0;
    }
    
    $log->info("Completed BOLD record updates: $updated_count records linked to taxonomy");
}

sub find_taxonid_for_record {
    my $record = shift;
    my $kingdom = $record->get_column('kingdom');
    
    # Build taxonomic path for this record
    my @path_components;
    for my $i (0 .. $#levels) {
        my $level = $levels[$i];
        my $name = $record->get_column($level) || '';
        
        # FIXED: Skip "None" values instead of stopping
        if ($name eq '' || $name eq 'None' || !defined($name)) {
            next;  # Skip this level but continue processing
        }
        
        push @path_components, {
            kingdom => $kingdom,
            name => $name,
            level => $level
        };
    }
    
    return undef unless @path_components;
    
    # Find taxonid for the most specific (last) taxon in the path
    my $lowest_taxon = $path_components[-1];
    
    # We need to find this taxon by walking down the taxonomic hierarchy
    my $parent_taxonid = undef;
    my $current_taxonid = undef;
    
    for my $component (@path_components) {
        my $sth = $dbh->prepare(
            "SELECT taxonid FROM taxa WHERE kingdom = ? AND name = ? AND level = ? AND " .
            "(parent_taxonid = ? OR (parent_taxonid IS NULL AND ? IS NULL))"
        );
        
        $sth->execute(
            $component->{kingdom},
            $component->{name}, 
            $component->{level},
            $parent_taxonid,
            $parent_taxonid
        );
        
        ($current_taxonid) = $sth->fetchrow_array;
        last unless $current_taxonid;
        
        $parent_taxonid = $current_taxonid;
    }
    
    return $current_taxonid;
}
