-- Calculate and store ranks and sumscores for all BOLD records
-- This script extracts the ranking logic from ranking_with_sumscore.sql 
-- and stores the results in the bold_ranks table for efficient querying

-- Clear any existing rank data
DELETE FROM bold_ranks;

-- Insert calculated ranks and sumscores for all records
INSERT INTO bold_ranks (recordid, ranking, sumscore)
WITH PivotedCriteria AS (
    SELECT
        bc.recordid,
        SUM(bc.status) AS sumscore,
        MAX(CASE WHEN c.name = 'SPECIES_ID' THEN bc.status ELSE NULL END) AS SPECIES_ID,
        MAX(CASE WHEN c.name = 'TYPE_SPECIMEN' THEN bc.status ELSE NULL END) AS TYPE_SPECIMEN,
        MAX(CASE WHEN c.name = 'SEQ_QUALITY' THEN bc.status ELSE NULL END) AS SEQ_QUALITY,
        MAX(CASE WHEN c.name = 'HAS_IMAGE' THEN bc.status ELSE NULL END) AS HAS_IMAGE,
        MAX(CASE WHEN c.name = 'COLLECTORS' THEN bc.status ELSE NULL END) AS COLLECTORS,
        MAX(CASE WHEN c.name = 'COLLECTION_DATE' THEN bc.status ELSE NULL END) AS COLLECTION_DATE,
        MAX(CASE WHEN c.name = 'COUNTRY' THEN bc.status ELSE NULL END) AS COUNTRY,
        MAX(CASE WHEN c.name = 'REGION' THEN bc.status ELSE NULL END) AS REGION,
        MAX(CASE WHEN c.name = 'SECTOR' THEN bc.status ELSE NULL END) AS SECTOR,
        MAX(CASE WHEN c.name = 'SITE' THEN bc.status ELSE NULL END) AS SITE,
        MAX(CASE WHEN c.name = 'COORD' THEN bc.status ELSE NULL END) AS COORD,
        MAX(CASE WHEN c.name = 'IDENTIFIER' THEN bc.status ELSE NULL END) AS IDENTIFIER,
        MAX(CASE WHEN c.name = 'ID_METHOD' THEN bc.status ELSE NULL END) AS ID_METHOD,
        MAX(CASE WHEN c.name = 'INSTITUTION' THEN bc.status ELSE NULL END) AS INSTITUTION,
        MAX(CASE WHEN c.name = 'PUBLIC_VOUCHER' THEN bc.status ELSE NULL END) AS PUBLIC_VOUCHER,
        MAX(CASE WHEN c.name = 'MUSEUM_ID' THEN bc.status ELSE NULL END) AS MUSEUM_ID,
        MAX(CASE WHEN c.name = 'HAPLOTYPE_ID' THEN bc.status ELSE NULL END) AS HAPLOTYPE_ID
    FROM
        bold_criteria bc
    JOIN
        criteria c ON bc.criterionid = c.criterionid
    GROUP BY
        bc.recordid
)
SELECT
    pc.recordid,
    CASE
        WHEN pc.SPECIES_ID = 1 AND pc.TYPE_SPECIMEN = 1 THEN 1
        WHEN pc.SPECIES_ID = 1 AND pc.SEQ_QUALITY = 1 AND pc.HAS_IMAGE = 1 AND pc.COLLECTORS = 1 AND pc.COLLECTION_DATE = 1 AND pc.COUNTRY = 1 AND (pc.SITE = 1 OR pc.SECTOR = 1 OR pc.REGION = 1 OR pc.COORD = 1) AND pc.IDENTIFIER = 1 AND pc.ID_METHOD = 1 AND (pc.INSTITUTION = 1 OR pc.PUBLIC_VOUCHER = 1 OR pc.MUSEUM_ID = 1) THEN 2
        WHEN pc.SPECIES_ID = 1 AND pc.SEQ_QUALITY = 1 AND pc.COUNTRY = 1 AND pc.IDENTIFIER = 1 AND pc.ID_METHOD = 1 AND (pc.INSTITUTION = 1 OR pc.PUBLIC_VOUCHER = 1 OR pc.MUSEUM_ID = 1) THEN 3
        WHEN pc.SPECIES_ID = 1 AND pc.SEQ_QUALITY = 1 AND pc.COUNTRY = 1 THEN 4
        WHEN pc.SPECIES_ID = 1 AND pc.SEQ_QUALITY = 1 THEN 5
        WHEN pc.SPECIES_ID = 1 THEN 6
        ELSE 7
    END AS ranking,
    pc.sumscore
FROM
    PivotedCriteria pc;

-- Create summary statistics for verification
SELECT 'Rank calculation completed' as status;
SELECT 'Total records processed: ' || COUNT(*) as summary FROM bold_ranks;
SELECT 'Rank distribution:' as summary;
SELECT 
    'Rank ' || ranking || ': ' || COUNT(*) || ' records' as distribution 
FROM bold_ranks 
GROUP BY ranking 
ORDER BY ranking;

SELECT 'Sumscore statistics:' as summary;
SELECT 
    'Min sumscore: ' || MIN(sumscore) || ', Max sumscore: ' || MAX(sumscore) || ', Avg sumscore: ' || ROUND(AVG(sumscore), 2) as stats
FROM bold_ranks;
