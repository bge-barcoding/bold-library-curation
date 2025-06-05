def create_single_family_db(family_info, source_db, output_file, logger):
    """Create database for a single family using the complete processed result data"""
    try:
        # Connect to source database and create target
        source_conn = sqlite3.connect(source_db)
        
        # Get the complete query from the ranking script
        ranking_query = """
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
        SELECT
            mc.url,
            mc.status,
            mc.additionalStatus,
            mc.curator_notes,
            br.sumscore,
            bags.bags_grade AS BAGS,
            br.ranking,
            CASE 
                WHEN cr.recordid IS NOT NULL THEN 'Yes'
                ELSE 'No'
            END AS country_representative,
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
            bo.otu_id
        FROM
            bold b
        LEFT JOIN
            manual_curation mc ON b.recordid = mc.recordid
        LEFT JOIN
            PivotedCriteria pc ON b.recordid = pc.recordid
        LEFT JOIN
            bags ON b.taxonid = bags.taxonid
        LEFT JOIN
            bold_haplotypes bh ON b.recordid = bh.recordid
        LEFT JOIN
            bold_otus bo ON b.recordid = bo.recordid
        LEFT JOIN
            bold_ranks br ON b.recordid = br.recordid
        LEFT JOIN
            country_representatives cr ON b.recordid = cr.recordid
        WHERE b.kingdom = ? AND b.phylum = ? AND b."class" = ? AND b."order" = ? AND b.family = ?
        """
        
        # Create new database with complete results as single table
        target_conn = sqlite3.connect(output_file)
        target_cursor = target_conn.cursor()
        
        # Execute query and create records table
        source_cursor = source_conn.cursor()
        source_cursor.execute(ranking_query, (family_info['kingdom'], family_info['phylum'], 
                                            family_info['class'], family_info['order'], family_info['family']))
        
        # Get column names
        columns = [description[0] for description in source_cursor.description]
        
        # Create table structure
        create_table_sql = f"CREATE TABLE records ({', '.join([f'{col} TEXT' for col in columns])})"
        target_cursor.execute(create_table_sql)
        
        # Insert all rows
        rows = source_cursor.fetchall()
        if rows:
            placeholders = ', '.join(['?' for _ in columns])
            target_cursor.executemany(f"INSERT INTO records VALUES ({placeholders})", rows)
            target_conn.commit()
            logger.info(f"Created family database: {output_file} with {len(rows)} records")
        else:
            logger.warning(f"No records found for family {family_info['family']}")
        
        source_conn.close()
        target_conn.close()
        return True
        
    except Exception as e:
        logger.error(f"Failed to create family database {output_file}: {e}")
        return False

def create_subfamily_db(family_info, subfamily, source_db, output_file, logger):
    """Create database for a specific subfamily using complete processed result data"""
    try:
        # Connect to source database
        source_conn = sqlite3.connect(source_db)
        
        # Use the same complete query but with subfamily filter
        ranking_query = """
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
        SELECT
            mc.url,
            mc.status,
            mc.additionalStatus,
            mc.curator_notes,
            br.sumscore,
            bags.bags_grade AS BAGS,
            br.ranking,
            CASE 
                WHEN cr.recordid IS NOT NULL THEN 'Yes'
                ELSE 'No'
            END AS country_representative,
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
            bo.otu_id
        FROM
            bold b
        LEFT JOIN
            manual_curation mc ON b.recordid = mc.recordid
        LEFT JOIN
            PivotedCriteria pc ON b.recordid = pc.recordid
        LEFT JOIN
            bags ON b.taxonid = bags.taxonid
        LEFT JOIN
            bold_haplotypes bh ON b.recordid = bh.recordid
        LEFT JOIN
            bold_otus bo ON b.recordid = bo.recordid
        LEFT JOIN
            bold_ranks br ON b.recordid = br.recordid
        LEFT JOIN
            country_representatives cr ON b.recordid = cr.recordid
        WHERE b.kingdom = ? AND b.phylum = ? AND b."class" = ? AND b."order" = ? AND b.family = ?
        """
        
        # Add subfamily condition
        if subfamily is None:
            ranking_query += " AND (b.subfamily IS NULL OR b.subfamily = '')"
            params = (family_info['kingdom'], family_info['phylum'], family_info['class'], 
                     family_info['order'], family_info['family'])
        else:
            ranking_query += " AND b.subfamily = ?"
            params = (family_info['kingdom'], family_info['phylum'], family_info['class'], 
                     family_info['order'], family_info['family'], subfamily)
        
        # Create new database
        target_conn = sqlite3.connect(output_file)
        target_cursor = target_conn.cursor()
        
        # Execute query
        source_cursor = source_conn.cursor()
        source_cursor.execute(ranking_query, params)
        
        # Get column names and create table
        columns = [description[0] for description in source_cursor.description]
        create_table_sql = f"CREATE TABLE records ({', '.join([f'{col} TEXT' for col in columns])})"
        target_cursor.execute(create_table_sql)
        
        # Insert data
        rows = source_cursor.fetchall()
        if rows:
            placeholders = ', '.join(['?' for _ in columns])
            target_cursor.executemany(f"INSERT INTO records VALUES ({placeholders})", rows)
            target_conn.commit()
            
            subfamily_label = subfamily if subfamily else "no_subfamily"
            logger.info(f"Created subfamily database: {output_file} ({subfamily_label}) with {len(rows)} records")
        else:
            subfamily_label = subfamily if subfamily else "no_subfamily"
            logger.warning(f"No records found for subfamily {subfamily_label}")
        
        source_conn.close()
        target_conn.close()
        return True
        
    except Exception as e:
        logger.error(f"Failed to create subfamily database {output_file}: {e}")
        return False
