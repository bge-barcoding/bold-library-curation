-- Select country representatives: Best record per species per OTU per country
-- Selection criteria: Lowest ranking, then highest sumscore, then lowest recordid (for deterministic tiebreaking)
-- Groups by: country_iso + species + otu_id

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

-- Generate summary statistics
SELECT 'Country representative selection completed' as status;

SELECT 'Summary Statistics:' as result;
SELECT 'Total country representatives selected: ' || COUNT(*) as summary FROM country_representatives;
SELECT 'Unique countries: ' || COUNT(DISTINCT country_iso) as summary FROM country_representatives;
SELECT 'Unique species: ' || COUNT(DISTINCT species) as summary FROM country_representatives;
SELECT 'Unique OTUs: ' || COUNT(DISTINCT otu_id) as summary FROM country_representatives;

SELECT 'Selection by ranking:' as result;
SELECT 
    'Rank ' || ranking || ': ' || COUNT(*) || ' representatives' as distribution 
FROM country_representatives 
GROUP BY ranking 
ORDER BY ranking;

SELECT 'Top 10 countries by number of representatives:' as result;
SELECT 
    country_iso, 
    COUNT(*) as representative_count,
    COUNT(DISTINCT species) as unique_species,
    COUNT(DISTINCT otu_id) as unique_otus
FROM country_representatives 
GROUP BY country_iso 
ORDER BY representative_count DESC 
LIMIT 10;

SELECT 'Top 10 species by number of country representatives:' as result;
SELECT 
    species, 
    COUNT(*) as country_count,
    COUNT(DISTINCT otu_id) as otu_count
FROM country_representatives 
GROUP BY species 
ORDER BY country_count DESC 
LIMIT 10;
