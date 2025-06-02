-- Select country representatives: Best record per species per OTU per country
-- Selection criteria: Lowest ranking, then highest sumscore, then lowest recordid (for deterministic tiebreaking)
-- Groups by: country_iso + species + otu_id
-- UPDATED: Excludes records with rank 7 OR BAGS grade F

-- Clear any existing country representative data
DELETE FROM country_representatives;

-- Insert country representatives using window functions for efficient selection
INSERT INTO country_representatives (recordid, country_iso, species, otu_id, ranking, sumscore, selection_reason)
WITH RankedRecords AS (
    SELECT 
        b.recordid,
        b.country_iso,
        b.species,
        bo.otu_id,
        br.ranking,
        br.sumscore,
        bags.bags_grade,
        -- Window function to rank records within each group
        ROW_NUMBER() OVER (
            PARTITION BY b.country_iso, b.species, bo.otu_id 
            ORDER BY br.ranking ASC, br.sumscore DESC, b.recordid ASC
        ) as selection_rank,
        -- Add selection reasoning
        CASE 
            WHEN ROW_NUMBER() OVER (
                PARTITION BY b.country_iso, b.species, bo.otu_id 
                ORDER BY br.ranking ASC, br.sumscore DESC, b.recordid ASC
            ) = 1 THEN 'Best quality record in group'
            ELSE 'Not selected'
        END as selection_reason
    FROM 
        bold b
    JOIN 
        bold_ranks br ON b.recordid = br.recordid
    JOIN 
        bold_otus bo ON b.recordid = bo.recordid
    LEFT JOIN 
        bags ON b.taxonid = bags.taxonid
    WHERE 
        -- Only include records with species-level identification
        b.species IS NOT NULL 
        AND b.species != ''
        AND b.species != 'None'          -- Exclude "None" values
        AND b.species NOT LIKE '%sp.%'   -- Exclude incomplete species identifications
        -- Only include records with country information
        AND b.country_iso IS NOT NULL 
        AND b.country_iso != ''
        AND b.country_iso != 'None'      -- Exclude "None" values
        -- Only include records assigned to OTUs
        AND bo.otu_id IS NOT NULL
        AND bo.otu_id != ''
        AND bo.otu_id != 'UNASSIGNED'
        -- NEW FILTERS: Exclude rank 7 OR BAGS grade F
        AND br.ranking != 7              -- Exclude lowest quality rank
        AND (bags.bags_grade IS NULL OR bags.bags_grade != 'F')  -- Exclude BAGS grade F (allow NULL for records without BAGS)
)
SELECT 
    recordid,
    country_iso,
    species,
    otu_id,
    ranking,
    sumscore,
    selection_reason
FROM RankedRecords 
WHERE selection_rank = 1;
