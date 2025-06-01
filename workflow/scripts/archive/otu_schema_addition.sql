-- Database schema updates to add OTU support
-- Add this to the existing schema_with_ranks.sql file

-- Table to store OTU assignments for BOLD records
-- Each record is assigned to exactly one OTU based on sequence similarity
CREATE TABLE IF NOT EXISTS "bold_otus" (
    "recordid" INTEGER NOT NULL,            -- Foreign key to bold table
    "otu_id" TEXT NOT NULL,                 -- OTU identifier (e.g., "OTU_000001")
    
    PRIMARY KEY (recordid),
    FOREIGN KEY(recordid) REFERENCES bold(recordid)
);

-- Indexes for OTU analysis performance
CREATE INDEX IF NOT EXISTS idx_bold_otus_recordid ON bold_otus(recordid);
CREATE INDEX IF NOT EXISTS idx_bold_otus_otu_id ON bold_otus(otu_id);
