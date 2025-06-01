-- Additional indexes for ranking and country representative selection performance
-- These indexes supplement the base indexes.sql file with rank-specific optimizations

-- Indexes for ranking performance
CREATE INDEX IF NOT EXISTS idx_bold_ranks_recordid ON bold_ranks(recordid);
CREATE INDEX IF NOT EXISTS idx_bold_ranks_ranking ON bold_ranks(ranking);
CREATE INDEX IF NOT EXISTS idx_bold_ranks_sumscore ON bold_ranks(sumscore);
CREATE INDEX IF NOT EXISTS idx_bold_ranks_ranking_sumscore ON bold_ranks(ranking, sumscore);

-- Indexes for country representative selection
CREATE INDEX IF NOT EXISTS idx_bold_country_species ON bold("country/ocean", species);
CREATE INDEX IF NOT EXISTS idx_bold_bin_uri ON bold(bin_uri);
CREATE INDEX IF NOT EXISTS idx_bold_species ON bold(species);
CREATE INDEX IF NOT EXISTS idx_bold_country_species_bin ON bold("country/ocean", species, bin_uri);

-- Composite index for efficient country representative queries
CREATE INDEX IF NOT EXISTS idx_bold_country_species_bin_haplotype ON bold("country/ocean", species, bin_uri) 
WHERE species IS NOT NULL AND species != '';

-- Index for haplotype-based grouping
CREATE INDEX IF NOT EXISTS idx_bold_haplotypes_haplotype_recordid ON bold_haplotypes(haplotype_id, recordid);
