-- Test script to verify subspecies BAGS inheritance
-- Run this after subspecies inheritance to check results

.headers ON
.mode column

-- Show subspecies that inherited BAGS grades
SELECT 
    'Subspecies with inherited BAGS grades:' as test_section;

SELECT 
    s.taxonid as subspecies_id,
    s.name as subspecies_name,
    parent.name as parent_species,
    b.bags_grade as inherited_grade,
    COUNT(b.bin_uri) as bin_count
FROM taxa s
JOIN taxa parent ON s.parent_taxonid = parent.taxonid
JOIN bags b ON s.taxonid = b.taxonid
WHERE s.level = 'subspecies' 
  AND parent.level = 'species'
GROUP BY s.taxonid, s.name, parent.name, b.bags_grade
ORDER BY b.bags_grade, parent.name, s.name
LIMIT 20;

-- Summary statistics
SELECT 
    'Summary statistics:' as test_section;

SELECT 
    'Total subspecies in taxa table:' as metric,
    COUNT(*) as count
FROM taxa 
WHERE level = 'subspecies';

SELECT 
    'Subspecies with inherited BAGS grades:' as metric,
    COUNT(DISTINCT s.taxonid) as count
FROM taxa s
JOIN bags b ON s.taxonid = b.taxonid
WHERE s.level = 'subspecies';

-- Grade distribution for subspecies
SELECT 
    'Grade distribution for subspecies:' as test_section;

SELECT 
    b.bags_grade,
    COUNT(*) as subspecies_count
FROM taxa s
JOIN bags b ON s.taxonid = b.taxonid
WHERE s.level = 'subspecies'
GROUP BY b.bags_grade
ORDER BY b.bags_grade;
