-- Additional table for storing country representatives
-- This table stores one representative record per species per OTU per country

CREATE TABLE IF NOT EXISTS "country_representatives" (
    "recordid" INTEGER NOT NULL,            -- Foreign key to bold table (the selected representative)
    "country_iso" TEXT NOT NULL,            -- Country ISO code
    "species" TEXT NOT NULL,                -- Species name
    "otu_id" TEXT NOT NULL,                 -- OTU identifier
    "ranking" INTEGER NOT NULL,             -- Ranking of the selected record
    "sumscore" INTEGER NOT NULL,            -- Sumscore of the selected record
    "selection_reason" TEXT,                -- Why this record was selected
    "selected_at" TEXT DEFAULT (datetime('now')), -- Timestamp when selection was made
    
    PRIMARY KEY (country_iso, species, otu_id),  -- One representative per country/species/OTU combination
    FOREIGN KEY(recordid) REFERENCES bold(recordid)
);

-- Indexes for country representative queries
CREATE INDEX IF NOT EXISTS idx_country_reps_recordid ON country_representatives(recordid);
CREATE INDEX IF NOT EXISTS idx_country_reps_country ON country_representatives(country_iso);
CREATE INDEX IF NOT EXISTS idx_country_reps_species ON country_representatives(species);
CREATE INDEX IF NOT EXISTS idx_country_reps_otu ON country_representatives(otu_id);
CREATE INDEX IF NOT EXISTS idx_country_reps_ranking ON country_representatives(ranking);
CREATE INDEX IF NOT EXISTS idx_country_reps_country_species ON country_representatives(country_iso, species);

-- Composite index for efficient lookups
CREATE INDEX IF NOT EXISTS idx_country_reps_composite ON country_representatives(country_iso, species, otu_id, ranking);
