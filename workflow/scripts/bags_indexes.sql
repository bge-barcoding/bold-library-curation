-- Additional indexes specifically for BAGS analysis optimization

-- Composite index for BIN sharing queries (most critical)
CREATE INDEX IF NOT EXISTS "bin_taxon_level_idx" ON bold ("bin_uri", "taxonid");

-- Index for subspecies name pattern matching
CREATE INDEX IF NOT EXISTS "taxa_level_name_idx" ON taxa ("level", "name");

-- Composite index for species-level taxa queries
CREATE INDEX IF NOT EXISTS "taxa_level_species_idx" ON taxa ("level") WHERE level = 'species';

-- Index for efficient BIN URI filtering
CREATE INDEX IF NOT EXISTS "bin_uri_valid_idx" ON bold ("bin_uri") WHERE bin_uri LIKE 'BOLD:%';

-- Index for joining bold and taxa tables efficiently
CREATE INDEX IF NOT EXISTS "bold_taxon_join_idx" ON bold ("taxonid", "bin_uri");

-- Index for taxonomic hierarchy traversal
CREATE INDEX IF NOT EXISTS "taxa_parent_level_idx" ON taxa ("parent_taxonid", "level");

-- Analyze tables to update query planner statistics
ANALYZE bold;
ANALYZE taxa;

-- Enable query planner optimizations
PRAGMA optimize;
PRAGMA journal_mode = WAL;  -- Enable Write-Ahead Logging for better concurrency
PRAGMA synchronous = NORMAL;  -- Balance safety vs performance
PRAGMA cache_size = -64000;   -- 64MB cache
PRAGMA temp_store = MEMORY;   -- Store temporary tables in memory
