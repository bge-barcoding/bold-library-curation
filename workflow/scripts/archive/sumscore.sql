-- This script operates on the database assuming all criteria have been assessed
-- and persisted into it. It calculates the sum of the status results for each
-- record, resulting in a `sumscore` column, whose values range from 0 to 16,
-- the latter being the number of pass/fail criteria that is assessed at time of
-- writing. Subsequently, the next query partitions the records into the constituent
-- bin_uri fields, sorts the records within each BIN by their sumscore in descending
-- order, emitting the highest scoring record for each BIN. To use this for producing
-- an output set in BCDM TSV format, do the following:
-- $ sqlite3 bold.db
-- sqlite> .headers ON
-- sqlite> .mode tabs
-- sqlite> .output curated.tsv
-- sqlite> .read sumscore.sql
-- sqlite> .quit

WITH ScoredRecords AS (
    SELECT
        b.*,
        SUM(bc.status) AS sumscore
    FROM
        bold b
    LEFT JOIN
        bold_criteria bc ON b.recordid = bc.recordid
    WHERE
        b.bin_uri <> ''
    GROUP BY
        b.recordid
)

SELECT
    sr.*
FROM (
    SELECT
        *,
        RANK() OVER (PARTITION BY bin_uri ORDER BY sumscore DESC) AS rank
    FROM
        ScoredRecords
) sr
WHERE
    sr.rank = 1;
