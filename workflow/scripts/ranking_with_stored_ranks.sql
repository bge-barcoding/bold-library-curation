-- Updated final output query that uses pre-calculated ranks from bold_ranks table
-- This replaces the complex CASE statement calculation with a simple JOIN
-- Provides better performance and consistency for large datasets

-- CTE to pivot the criteria into separate columns based on recordid
WITH PivotedCriteria AS (
    SELECT
        bc.recordid,
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
-- Main query to select from bold and join on pivoted criteria, BAGS grades, and pre-calculated ranks
SELECT
    b.*,
    pc.SPECIES_ID,
    pc.TYPE_SPECIMEN,
    pc.SEQ_QUALITY,
    pc.HAS_IMAGE,
    pc.COLLECTORS AS HAS_COLLECTORS,
    pc.COLLECTION_DATE,
    pc.COUNTRY AS HAS_COUNTRY,
    pc.REGION AS HAS_REGION,
    pc.SECTOR AS HAS_SECTOR,
    pc.SITE AS HAS_SITE,
    pc.COORD AS HAS_COORD,
    pc.IDENTIFIER,
    pc.ID_METHOD,
    pc.INSTITUTION,
    pc.PUBLIC_VOUCHER,
    pc.MUSEUM_ID,
    bh.haplotype_id,
    br.sumscore,        -- Use pre-calculated sumscore from bold_ranks
    bags.bags_grade,
    br.ranking          -- Use pre-calculated ranking from bold_ranks
FROM
    bold b
LEFT JOIN
    PivotedCriteria pc ON b.recordid = pc.recordid
LEFT JOIN
    bags ON b.taxonid = bags.taxonid
LEFT JOIN
    bold_haplotypes bh ON b.recordid = bh.recordid
LEFT JOIN
    bold_ranks br ON b.recordid = br.recordid;  -- Join with pre-calculated ranks
