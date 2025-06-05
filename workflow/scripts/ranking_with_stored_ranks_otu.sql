-- Updated final output query that uses pre-calculated ranks from bold_ranks table
-- This replaces the complex CASE statement calculation with a simple JOIN
-- Provides better performance and consistency for large datasets
-- UPDATED: Includes OTU_ID from bold_otus table

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
-- Main query to select from bold and join on pivoted criteria, BAGS grades, haplotypes, OTUs, country representatives, pre-calculated ranks, and manual curation
SELECT
    -- Manual curation fields at the front
    mc.url,
    mc.status,
    mc.additionalStatus,
    mc.curator_notes,
    -- Key ranking fields
    br.sumscore,            -- Use pre-calculated sumscore from bold_ranks
    bags.bags_grade AS bags_static,
	CASE 
		WHEN bags.bags_grade = 'C' THEN 'BIN splitting: ' || bags.bin_uri
		WHEN bags.bags_grade = 'E' THEN 'BIN sharing: ' || bags.bin_uri || '; ' || bags.sharers
		ELSE 'Single BIN: ' || bags.bin_uri
	END AS "BIN info",
    br.ranking,             -- Use pre-calculated ranking from bold_ranks
    CASE 
        WHEN cr.recordid IS NOT NULL THEN 'Yes'
        ELSE 'No'
    END AS country_representative,  -- Indicate if record is a country representative
    -- Original bold table fields (with renamed columns)
    b.recordid,
    b.taxonid,
    b.processid,
    b.sampleid,
    b.fieldid,
    b.museumid,
    b.record_id,
    b.specimenid,
    b.processid_minted_date,
    b.bin_uri,
    b.bin_created_date,
    b.collection_code,
    b.inst,
    b.taxid,
    b.taxon_name,
    b.taxon_rank,
    b.kingdom,
    b.phylum,
    b."class",
    b."order",
    b.family,
    b.subfamily,
    b.tribe,
    b.genus,
    b.species,
    b.subspecies,
    b.species_reference,
    b.identification,
    b.identification_method,
    b.identification_rank,
    b.identified_by,
    b.identifier_email,
    b.taxonomy_notes,
    b.sex,
    b.reproduction,
    b.life_stage,
    b.short_note,
    b.notes,
    b.voucher_type,
    b.tissue_type,
    b.specimen_linkout,
    b.associated_specimens,
    b.associated_taxa,
    b.collectors,
    b.collection_date_start,
    b.collection_date_end,
    b.collection_event_id,
    b.collection_time,
    b.collection_notes,
    b.geoid,
    b."country/ocean" AS country_ocean,
    b.country_iso,
    b."province/state" AS province_state,
    b.region,
    b.sector,
    b.site,
    b.site_code,
    b.coord,
    b.coord_accuracy,
    b.coord_source,
    b.elev,
    b.elev_accuracy,
    b."depth",
    b.depth_accuracy,
    b.habitat,
    b.sampling_protocol,
    b.nuc,
    b.nuc_basecount,
    b.insdc_acs,
    b.funding_src,
    b.marker_code,
    b.primers_forward,
    b.primers_reverse,
    b.sequence_run_site,
    b.sequence_upload_date,
    b.bold_recordset_code_arr AS recordset_code_arr,
    b.sovereign_inst,
    b.realm,
    b.biome,
    b.ecoregion,
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
    bo.otu_id              -- Include OTU_ID from bold_otus table
FROM
    bold b
LEFT JOIN
    manual_curation mc ON b.recordid = mc.recordid    -- Join with manual curation data
LEFT JOIN
    PivotedCriteria pc ON b.recordid = pc.recordid
LEFT JOIN
    bags ON b.taxonid = bags.taxonid
LEFT JOIN
    bold_haplotypes bh ON b.recordid = bh.recordid
LEFT JOIN
    bold_otus bo ON b.recordid = bo.recordid          -- Join with OTU assignments
LEFT JOIN
    bold_ranks br ON b.recordid = br.recordid         -- Join with pre-calculated ranks
LEFT JOIN
    country_representatives cr ON b.recordid = cr.recordid;  -- NEW: Join with country representatives
