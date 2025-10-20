-- Critical composite indexes for country representative selection
-- These must be created BEFORE running select_country_representatives.sql

-- 1. Main composite index for window function partitioning
CREATE INDEX IF NOT EXISTS idx_bold_country_species_otu 
ON bold(country_iso, species);

-- 2. Composite index for OTU joins with filtering
CREATE INDEX IF NOT EXISTS idx_bold_otus_recordid_otu 
ON bold_otus(recordid, otu_id) 
WHERE otu_id IS NOT NULL 
AND otu_id != '' 
AND otu_id != 'UNASSIGNED';

-- 3. Critical index for EXISTS subquery on bold_criteria
CREATE INDEX IF NOT EXISTS idx_bold_criteria_species_check 
ON bold_criteria(recordid, criterionid, status) 
WHERE criterionid = 1 AND status = 1;

-- 4. Composite index for bold filtering
CREATE INDEX IF NOT EXISTS idx_bold_filters 
ON bold(recordid, country_iso, bin_uri, species, taxonid) 
WHERE country_iso IS NOT NULL 
AND country_iso != '' 
AND country_iso != 'None'
AND bin_uri IS NOT NULL 
AND bin_uri != 'None';

-- 5. Index for bags join with grade F filtering
CREATE INDEX IF NOT EXISTS idx_bags_taxonid_grade 
ON bags(taxonid, bags_grade);

-- 6. Composite index for bold_ranks filtering
CREATE INDEX IF NOT EXISTS idx_bold_ranks_filter 
ON bold_ranks(recordid, ranking) 
WHERE ranking != 7;

-- Analyze tables after creating indexes
ANALYZE bold;
ANALYZE bold_otus;
ANALYZE bold_criteria;
ANALYZE bold_ranks;
ANALYZE bags;

-- Optimize query planner
PRAGMA optimize;