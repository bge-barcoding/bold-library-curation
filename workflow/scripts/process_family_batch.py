#!/usr/bin/env python3
"""
Family Batch Processor - Processes a batch of families in parallel
Designed for SLURM job array execution
"""

import sqlite3
import argparse
import json
import os
import sys
import multiprocessing as mp
from pathlib import Path
import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
import traceback
import time

def get_complete_ranking_query():
    """Return the complete ranking query from ranking_with_stored_ranks_otu.sql"""
    return """
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
        CASE 
            WHEN bags.bags_grade = 'C' THEN 'BIN splitting: ' || bags.bin_uri
            WHEN bags.bags_grade = 'E' THEN 'BIN sharing: ' || bags.bin_uri || '; ' || bags.sharers
            ELSE 'Single BIN: ' || bags.bin_uri
        END AS "BIN info",
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
        country_representatives cr ON b.recordid = cr.recordid"""

def setup_logging(log_file):
    """Setup logging for the batch processor"""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def create_family_database(family_info, source_db, output_dir, threshold, logger, retry_count=0):
    """Create database for a single family with subfamily splitting if needed"""
    try:
        # Build taxonomy path
        kingdom = family_info['kingdom'].replace('/', '_').replace(' ', '_')
        phylum = family_info['phylum'].replace('/', '_').replace(' ', '_')
        class_name = family_info['class'].replace('/', '_').replace(' ', '_')
        order = family_info['order'].replace('/', '_').replace(' ', '_')
        family = family_info['family'].replace('/', '_').replace(' ', '_')
        
        # Create directory structure
        family_dir = Path(output_dir) / kingdom / phylum / class_name / order / family
        family_dir.mkdir(parents=True, exist_ok=True)
        
        conn = sqlite3.connect(source_db)
        cursor = conn.cursor()
        
        # Check if family needs splitting by subfamily
        if family_info['needs_splitting']:
            logger.info(f"Processing large family {family} with subfamily splitting")
            
            # Get subfamilies for this family
            subfamily_query = """
            SELECT DISTINCT subfamily, COUNT(*) as count
            FROM bold 
            WHERE kingdom = ? AND phylum = ? AND "class" = ? AND "order" = ? AND family = ?
            AND subfamily IS NOT NULL AND subfamily != ''
            GROUP BY subfamily
            """
            cursor.execute(subfamily_query, (family_info['kingdom'], family_info['phylum'], 
                                           family_info['class'], family_info['order'], family_info['family']))
            subfamilies = cursor.fetchall()
            
            if subfamilies:
                created_files = []
                for subfamily, count in subfamilies:
                    subfamily_clean = subfamily.replace('/', '_').replace(' ', '_')
                    db_file = family_dir / f"{family}_{subfamily_clean}.db"
                    
                    if create_subfamily_db(family_info, subfamily, source_db, db_file, logger):
                        created_files.append(str(db_file))
                
                # Handle records without subfamily
                no_subfamily_query = """
                SELECT COUNT(*) FROM bold 
                WHERE kingdom = ? AND phylum = ? AND "class" = ? AND "order" = ? AND family = ?
                AND (subfamily IS NULL OR subfamily = '')
                """
                cursor.execute(no_subfamily_query, (family_info['kingdom'], family_info['phylum'], 
                                                  family_info['class'], family_info['order'], family_info['family']))
                no_subfamily_count = cursor.fetchone()[0]
                
                if no_subfamily_count > 0:
                    db_file = family_dir / f"{family}.db"
                    if create_subfamily_db(family_info, None, source_db, db_file, logger):
                        created_files.append(str(db_file))
                
                conn.close()
                return created_files
        
        # Create single family database
        db_file = family_dir / f"{family}.db"
        success = create_single_family_db(family_info, source_db, db_file, logger)
        conn.close()
        
        return [str(db_file)] if success else []
        
    except Exception as e:
        if retry_count < 2:  # Retry up to 2 times
            logger.warning(f"Retrying family {family_info['family']} (attempt {retry_count + 1}): {e}")
            time.sleep(1)  # Brief delay before retry
            return create_family_database(family_info, source_db, output_dir, threshold, logger, retry_count + 1)
        else:
            logger.error(f"Failed to process family {family_info['family']} after retries: {e}")
            logger.error(traceback.format_exc())
            return []

def create_single_family_db(family_info, source_db, output_file, logger):
    """Create database for a single family using complete processed result data"""
    try:
        # Connect to source database
        source_conn = sqlite3.connect(source_db)
        
        # Use the complete ranking query (same as ranking_with_stored_ranks_otu.sql)
        ranking_query = get_complete_ranking_query()
        ranking_query += " WHERE b.kingdom = ? AND b.phylum = ? AND b.\"class\" = ? AND b.\"order\" = ? AND b.family = ?"
        
        # Create new database with complete results as single table
        target_conn = sqlite3.connect(output_file)
        target_cursor = target_conn.cursor()
        
        # Execute query
        source_cursor = source_conn.cursor()
        source_cursor.execute(ranking_query, (family_info['kingdom'], family_info['phylum'], 
                                            family_info['class'], family_info['order'], family_info['family']))
        
        # Get column names and create table
        columns = [description[0] for description in source_cursor.description]
        create_table_sql = f"CREATE TABLE records ({', '.join([f'[{col}] TEXT' for col in columns])})"
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
        source_conn = sqlite3.connect(source_db)
        
        # Use the complete ranking query with subfamily filter
        ranking_query = get_complete_ranking_query()
        
        if subfamily is None:
            ranking_query += " WHERE b.kingdom = ? AND b.phylum = ? AND b.\"class\" = ? AND b.\"order\" = ? AND b.family = ? AND (b.subfamily IS NULL OR b.subfamily = '')"
            params = (family_info['kingdom'], family_info['phylum'], family_info['class'], 
                     family_info['order'], family_info['family'])
        else:
            ranking_query += " WHERE b.kingdom = ? AND b.phylum = ? AND b.\"class\" = ? AND b.\"order\" = ? AND b.family = ? AND b.subfamily = ?"
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
        create_table_sql = f"CREATE TABLE records ({', '.join([f'[{col}] TEXT' for col in columns])})"
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


def process_batch(batch_file, source_db, output_dir, threshold, max_workers, logger):
    """Process a batch of families using thread pool"""
    try:
        with open(batch_file, 'r') as f:
            families = json.load(f)
        
        if not families:
            logger.info("Empty batch - nothing to process")
            return {'processed': 0, 'failed': 0, 'created_files': []}
        
        logger.info(f"Processing batch with {len(families)} families using {max_workers} workers")
        
        processed = 0
        failed = 0
        created_files = []
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Submit all family processing tasks
            future_to_family = {
                executor.submit(create_family_database, family, source_db, output_dir, threshold, logger): family
                for family in families
            }
            
            # Process completed tasks
            for future in as_completed(future_to_family):
                family = future_to_family[future]
                try:
                    result = future.result()
                    if result:
                        processed += 1
                        created_files.extend(result)
                        logger.info(f"✓ Completed family: {family['family']} ({len(result)} files)")
                    else:
                        failed += 1
                        logger.error(f"✗ Failed family: {family['family']}")
                except Exception as e:
                    failed += 1
                    logger.error(f"✗ Exception processing family {family['family']}: {e}")
        
        logger.info(f"Batch completed: {processed} processed, {failed} failed")
        return {
            'processed': processed,
            'failed': failed,
            'created_files': created_files
        }
        
    except Exception as e:
        logger.error(f"Error processing batch {batch_file}: {e}")
        logger.error(traceback.format_exc())
        return {'processed': 0, 'failed': 0, 'created_files': []}

def main():
    parser = argparse.ArgumentParser(description='Process a batch of families for database splitting')
    parser.add_argument('batch_file', help='JSON file containing family batch')
    parser.add_argument('--source-db', required=True, help='Source BOLD database')
    parser.add_argument('--output-dir', required=True, help='Output directory for family databases')
    parser.add_argument('--threshold', type=int, default=10000, help='Family size threshold')
    parser.add_argument('--max-workers', type=int, help='Maximum worker threads (default: CPU count)')
    parser.add_argument('--log-file', help='Log file path')
    
    args = parser.parse_args()
    
    # Set default workers
    if args.max_workers is None:
        args.max_workers = min(mp.cpu_count(), 8)  # Cap at 8 to avoid overwhelming
    
    # Setup logging
    if args.log_file:
        logger = setup_logging(args.log_file)
    else:
        batch_name = Path(args.batch_file).stem
        log_file = Path(args.output_dir) / f"{batch_name}.log"
        logger = setup_logging(log_file)
    
    logger.info(f"Starting batch processor")
    logger.info(f"Batch file: {args.batch_file}")
    logger.info(f"Source database: {args.source_db}")
    logger.info(f"Output directory: {args.output_dir}")
    logger.info(f"Max workers: {args.max_workers}")
    
    # Process the batch
    start_time = time.time()
    result = process_batch(args.batch_file, args.source_db, args.output_dir, args.threshold, args.max_workers, logger)
    end_time = time.time()
    
    # Write results summary
    batch_name = Path(args.batch_file).stem
    result_file = Path(args.output_dir) / f"{batch_name}_result.json"
    
    result['batch_file'] = args.batch_file
    result['processing_time'] = end_time - start_time
    result['timestamp'] = time.strftime('%Y-%m-%d %H:%M:%S')
    
    with open(result_file, 'w') as f:
        json.dump(result, f, indent=2)
    
    logger.info(f"Results written to: {result_file}")
    logger.info(f"Processing time: {result['processing_time']:.2f} seconds")
    
    # Exit with error code if there were failures
    if result['failed'] > 0:
        logger.warning(f"Completed with {result['failed']} failures")
        sys.exit(1)
    else:
        logger.info("All families processed successfully")

if __name__ == "__main__":
    main()
