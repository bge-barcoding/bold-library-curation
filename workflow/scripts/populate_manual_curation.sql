-- Populate manual_curation table with URLs for all BOLD records
-- This script should run after all assessments but before TSV export
-- Handles cases where processid might be NULL or empty

INSERT INTO manual_curation (recordid, processid, url)
SELECT 
    recordid, 
    processid,
    CASE 
        WHEN processid IS NULL OR processid = '' THEN 'https://portal.boldsystems.org/record/None'
        ELSE 'https://portal.boldsystems.org/record/' || processid
    END as url
FROM bold 
WHERE recordid NOT IN (SELECT recordid FROM manual_curation);
