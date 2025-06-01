-- To execute me on the database to produce a TSV file, do:
-- $ sqlite3 mydatabase.db
-- sqlite> .headers ON
-- sqlite> .mode tabs
-- sqlite> .output result_output.tsv
-- sqlite> .read pivot.sql
-- sqlite> .quit


SELECT
  b.*,
  MAX(CASE WHEN c.name = 'SPECIES_ID' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "SPECIES_ID",
  MAX(CASE WHEN c.name = 'TYPE_SPECIMEN' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "TYPE_SPECIMEN",
  MAX(CASE WHEN c.name = 'SEQ_QUALITY' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "SEQ_QUALITY",
  MAX(CASE WHEN c.name = 'PUBLIC_VOUCHER' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "PUBLIC_VOUCHER",
  MAX(CASE WHEN c.name = 'HAS_IMAGE' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "HAS_IMAGE",
  MAX(CASE WHEN c.name = 'IDENTIFIER' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "IDENTIFIER",
  MAX(CASE WHEN c.name = 'ID_METHOD' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "ID_METHOD",
  MAX(CASE WHEN c.name = 'COLLECTORS' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "COLLECTORS",
  MAX(CASE WHEN c.name = 'COLLECTION_DATE' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "COLLECTION_DATE",
  MAX(CASE WHEN c.name = 'COUNTRY' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "COUNTRY",
  MAX(CASE WHEN c.name = 'SITE' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "SITE",
  MAX(CASE WHEN c.name = 'COORD' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "COORD",
  MAX(CASE WHEN c.name = 'INSTITUTION' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "INSTITUTION",
  MAX(CASE WHEN c.name = 'MUSEUM_ID' THEN CASE bc.status WHEN 1 THEN 'pass' ELSE 'fail' END END) AS "MUSEUM_ID"
FROM
  bold b
LEFT JOIN
  bold_criteria bc ON b.recordid = bc.recordid
LEFT JOIN
  criteria c ON bc.criterionid = c.criterionid
GROUP BY
  b.recordid;
