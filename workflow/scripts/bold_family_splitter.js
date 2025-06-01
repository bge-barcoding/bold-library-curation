#!/usr/bin/env node

/**
 * BOLD Family Database Splitter
 * 
 * Splits BOLD database or TSV export into family-level databases
 * organized by taxonomic hierarchy with configurable size thresholds
 */

const fs = require('fs');
const path = require('path');
const sqlite3 = require('sqlite3').verbose();
const { createReadStream } = require('fs');
const { createInterface } = require('readline');

// Configuration
const CONFIG = {
    // Input files (can be either database or TSV)
    INPUT_DB: null,
    INPUT_TSV: null,
    
    // Output configuration
    OUTPUT_BASE: 'taxonomic_output',
    
    // Thresholds
    FAMILY_SIZE_THRESHOLD: 10000,  // Split families larger than this by subfamily
    
    // Fallback values for missing taxonomy
    UNKNOWN_PHYLUM: 'Unknown_Phylum',
    UNKNOWN_CLASS: 'Unknown_Class', 
    UNKNOWN_ORDER: 'Unknown_Order',
    UNKNOWN_FAMILY: 'Unknown_Family',
    UNKNOWN_SUBFAMILY: 'Unknown_Subfamily'
};

class BoldFamilySplitter {
    constructor(inputPath, outputBase = CONFIG.OUTPUT_BASE) {
        this.inputPath = inputPath;
        this.outputBase = outputBase;
        this.isDatabase = inputPath.endsWith('.db');
        this.familyStats = new Map();
        this.taxonomyCache = new Map();
        
        // Ensure output directory exists
        fs.mkdirSync(this.outputBase, { recursive: true });
    }

    async run() {
        console.log(`Processing ${this.isDatabase ? 'database' : 'TSV'}: ${this.inputPath}`);
        console.log(`Output directory: ${this.outputBase}`);
        console.log(`Family size threshold: ${CONFIG.FAMILY_SIZE_THRESHOLD} records`);
        
        try {
            // Step 1: Analyze family sizes and get taxonomy info
            await this.analyzeFamilies();
            
            // Step 2: Create family databases with proper folder structure
            await this.createFamilyDatabases();
            
            // Step 3: Generate summary report
            this.generateReport();
            
            console.log('\n✅ Processing completed successfully!');
            
        } catch (error) {
            console.error('❌ Error:', error.message);
            throw error;
        }
    }

    async analyzeFamilies() {
        console.log('\n📊 Analyzing family sizes and taxonomy...');
        
        if (this.isDatabase) {
            await this.analyzeFamiliesFromDB();
        } else {
            await this.analyzeFamiliesFromTSV();
        }
        
        console.log('\nFamily analysis results:');
        for (const [family, stats] of this.familyStats) {
            const strategy = stats.count > CONFIG.FAMILY_SIZE_THRESHOLD ? 'subfamily' : 'family';
            console.log(`  ${family}: ${stats.count} records (split by ${strategy})`);
        }
    }

    async analyzeFamiliesFromDB() {
        const db = new sqlite3.Database(this.inputPath, sqlite3.OPEN_READONLY);
        
        try {
            // Get family counts
            const familyCounts = await this.queryAsync(db, `
                SELECT 
                    COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') as family,
                    COUNT(*) as count
                FROM records 
                GROUP BY family
                ORDER BY count DESC
            `);
            
            // Get taxonomy information for each family
            for (const row of familyCounts) {
                const family = row.family;
                
                // Get sample taxonomy for this family
                const taxonomyRow = await this.queryAsync(db, `
                    SELECT 
                        COALESCE(phylum, '${CONFIG.UNKNOWN_PHYLUM}') as phylum,
                        COALESCE(class, '${CONFIG.UNKNOWN_CLASS}') as class,
                        COALESCE("order", '${CONFIG.UNKNOWN_ORDER}') as order_name,
                        COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') as family,
                        subfamily
                    FROM records 
                    WHERE COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') = ?
                    LIMIT 1
                `, [family]);
                
                const taxonomy = taxonomyRow[0] || {};
                
                // Get subfamily counts if family is large
                let subfamilies = [];
                if (row.count > CONFIG.FAMILY_SIZE_THRESHOLD) {
                    const subfamilyRows = await this.queryAsync(db, `
                        SELECT 
                            COALESCE(subfamily, '${CONFIG.UNKNOWN_SUBFAMILY}') as subfamily,
                            COUNT(*) as count
                        FROM records 
                        WHERE COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') = ?
                        GROUP BY subfamily
                        ORDER BY count DESC
                    `, [family]);
                    
                    subfamilies = subfamilyRows;
                }
                
                this.familyStats.set(family, {
                    count: row.count,
                    taxonomy: {
                        phylum: taxonomy.phylum || CONFIG.UNKNOWN_PHYLUM,
                        class: taxonomy.class || CONFIG.UNKNOWN_CLASS,
                        order: taxonomy.order_name || CONFIG.UNKNOWN_ORDER,
                        family: family
                    },
                    subfamilies: subfamilies
                });
            }
            
        } finally {
            db.close();
        }
    }

    async analyzeFamiliesFromTSV() {
        const familyCounts = new Map();
        const taxonomyExamples = new Map();
        
        const fileStream = createReadStream(this.inputPath);
        const rl = createInterface({ input: fileStream, crlfDelay: Infinity });
        
        let headers = [];
        let lineCount = 0;
        
        // Find column indices
        let familyCol, phylumCol, classCol, orderCol, subfamilyCol;
        
        for await (const line of rl) {
            if (lineCount === 0) {
                headers = line.split('\t');
                
                // Find relevant columns (case-insensitive)
                familyCol = this.findColumnIndex(headers, ['family']);
                phylumCol = this.findColumnIndex(headers, ['phylum']);
                classCol = this.findColumnIndex(headers, ['class']);
                orderCol = this.findColumnIndex(headers, ['order']);
                subfamilyCol = this.findColumnIndex(headers, ['subfamily']);
                
                if (familyCol === -1) {
                    throw new Error('Family column not found in TSV');
                }
                
                lineCount++;
                continue;
            }
            
            const columns = line.split('\t');
            const family = columns[familyCol] || CONFIG.UNKNOWN_FAMILY;
            
            // Count families
            familyCounts.set(family, (familyCounts.get(family) || 0) + 1);
            
            // Store taxonomy example for each family
            if (!taxonomyExamples.has(family)) {
                taxonomyExamples.set(family, {
                    phylum: (phylumCol !== -1 ? columns[phylumCol] : null) || CONFIG.UNKNOWN_PHYLUM,
                    class: (classCol !== -1 ? columns[classCol] : null) || CONFIG.UNKNOWN_CLASS,
                    order: (orderCol !== -1 ? columns[orderCol] : null) || CONFIG.UNKNOWN_ORDER,
                    subfamily: subfamilyCol !== -1 ? columns[subfamilyCol] : null
                });
            }
            
            lineCount++;
            if (lineCount % 100000 === 0) {
                process.stdout.write(`\rProcessed ${lineCount} records...`);
            }
        }
        
        console.log(`\nProcessed ${lineCount - 1} total records`);
        
        // Analyze subfamily distribution for large families
        for (const [family, count] of familyCounts) {
            const taxonomy = taxonomyExamples.get(family);
            let subfamilies = [];
            
            if (count > CONFIG.FAMILY_SIZE_THRESHOLD) {
                // Re-scan to get subfamily counts for this family
                subfamilies = await this.getSubfamilyCountsFromTSV(family, familyCol, subfamilyCol);
            }
            
            this.familyStats.set(family, {
                count,
                taxonomy,
                subfamilies
            });
        }
    }

    async getSubfamilyCountsFromTSV(targetFamily, familyCol, subfamilyCol) {
        if (subfamilyCol === -1) return [];
        
        const subfamilyCounts = new Map();
        const fileStream = createReadStream(this.inputPath);
        const rl = createInterface({ input: fileStream, crlfDelay: Infinity });
        
        let lineCount = 0;
        
        for await (const line of rl) {
            if (lineCount === 0) {
                lineCount++;
                continue;
            }
            
            const columns = line.split('\t');
            const family = columns[familyCol] || CONFIG.UNKNOWN_FAMILY;
            
            if (family === targetFamily) {
                const subfamily = columns[subfamilyCol] || CONFIG.UNKNOWN_SUBFAMILY;
                subfamilyCounts.set(subfamily, (subfamilyCounts.get(subfamily) || 0) + 1);
            }
            
            lineCount++;
        }
        
        return Array.from(subfamilyCounts.entries()).map(([subfamily, count]) => ({
            subfamily,
            count
        }));
    }

    findColumnIndex(headers, possibleNames) {
        for (const name of possibleNames) {
            const index = headers.findIndex(h => h.toLowerCase() === name.toLowerCase());
            if (index !== -1) return index;
        }
        return -1;
    }

    async createFamilyDatabases() {
        console.log('\n🗂️  Creating family databases...');
        
        for (const [family, stats] of this.familyStats) {
            const { taxonomy, subfamilies } = stats;
            
            // Create taxonomic folder structure
            const familyDir = path.join(
                this.outputBase,
                this.sanitizeFileName(taxonomy.phylum),
                this.sanitizeFileName(taxonomy.class),
                this.sanitizeFileName(taxonomy.order)
            );
            
            fs.mkdirSync(familyDir, { recursive: true });
            
            if (subfamilies.length > 0) {
                // Split by subfamily
                console.log(`  Creating subfamily databases for ${family}...`);
                await this.createSubfamilyDatabases(family, familyDir, subfamilies);
            } else {
                // Create single family database
                const dbPath = path.join(familyDir, `${this.sanitizeFileName(family)}.db`);
                await this.createSingleFamilyDatabase(family, dbPath);
                console.log(`  ✅ Created ${dbPath}`);
            }
        }
    }

    async createSubfamilyDatabases(family, familyDir, subfamilies) {
        const subfamilyDir = path.join(familyDir, this.sanitizeFileName(family));
        fs.mkdirSync(subfamilyDir, { recursive: true });
        
        for (const subfamilyInfo of subfamilies) {
            const subfamily = subfamilyInfo.subfamily;
            const dbPath = path.join(subfamilyDir, `${this.sanitizeFileName(subfamily)}.db`);
            
            if (this.isDatabase) {
                await this.createSubfamilyDatabaseFromDB(family, subfamily, dbPath);
            } else {
                await this.createSubfamilyDatabaseFromTSV(family, subfamily, dbPath);
            }
            
            console.log(`    ✅ Created ${dbPath} (${subfamilyInfo.count} records)`);
        }
        
        // Create shared BINs analysis for this family
        await this.analyzeSharedBINs(family, subfamilyDir);
    }

    async createSingleFamilyDatabase(family, dbPath) {
        if (this.isDatabase) {
            await this.createFamilyDatabaseFromDB(family, dbPath);
        } else {
            await this.createFamilyDatabaseFromTSV(family, dbPath);
        }
    }

    async createFamilyDatabaseFromDB(family, dbPath) {
        const sourceDb = new sqlite3.Database(this.inputPath, sqlite3.OPEN_READONLY);
        const targetDb = new sqlite3.Database(dbPath);
        
        try {
            // Copy schema
            const schema = await this.queryAsync(sourceDb, "SELECT sql FROM sqlite_master WHERE type='table' AND name='records'");
            if (schema.length > 0) {
                await this.execAsync(targetDb, schema[0].sql);
            }
            
            // Copy data
            await this.execAsync(targetDb, 'ATTACH DATABASE ? AS source', [this.inputPath]);
            await this.execAsync(targetDb, `
                INSERT INTO records 
                SELECT * FROM source.records 
                WHERE COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') = ?
            `, [family]);
            
            // Create indexes
            await this.createIndexes(targetDb);
            
        } finally {
            sourceDb.close();
            targetDb.close();
        }
    }

    async createSubfamilyDatabaseFromDB(family, subfamily, dbPath) {
        const sourceDb = new sqlite3.Database(this.inputPath, sqlite3.OPEN_READONLY);
        const targetDb = new sqlite3.Database(dbPath);
        
        try {
            // Copy schema
            const schema = await this.queryAsync(sourceDb, "SELECT sql FROM sqlite_master WHERE type='table' AND name='records'");
            if (schema.length > 0) {
                await this.execAsync(targetDb, schema[0].sql);
            }
            
            // Copy data
            await this.execAsync(targetDb, 'ATTACH DATABASE ? AS source', [this.inputPath]);
            await this.execAsync(targetDb, `
                INSERT INTO records 
                SELECT * FROM source.records 
                WHERE COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') = ? 
                AND COALESCE(subfamily, '${CONFIG.UNKNOWN_SUBFAMILY}') = ?
            `, [family, subfamily]);
            
            // Create indexes
            await this.createIndexes(targetDb);
            
        } finally {
            sourceDb.close();
            targetDb.close();
        }
    }

    async createFamilyDatabaseFromTSV(family, dbPath) {
        const db = new sqlite3.Database(dbPath);
        
        try {
            // Create table structure
            await this.createTableFromTSV(db);
            
            // Insert data
            await this.insertTSVData(db, (columns, headers) => {
                const familyCol = this.findColumnIndex(headers, ['family']);
                const recordFamily = columns[familyCol] || CONFIG.UNKNOWN_FAMILY;
                return recordFamily === family;
            });
            
            // Create indexes
            await this.createIndexes(db);
            
        } finally {
            db.close();
        }
    }

    async createSubfamilyDatabaseFromTSV(family, subfamily, dbPath) {
        const db = new sqlite3.Database(dbPath);
        
        try {
            // Create table structure
            await this.createTableFromTSV(db);
            
            // Insert data
            await this.insertTSVData(db, (columns, headers) => {
                const familyCol = this.findColumnIndex(headers, ['family']);
                const subfamilyCol = this.findColumnIndex(headers, ['subfamily']);
                
                const recordFamily = columns[familyCol] || CONFIG.UNKNOWN_FAMILY;
                const recordSubfamily = columns[subfamilyCol] || CONFIG.UNKNOWN_SUBFAMILY;
                
                return recordFamily === family && recordSubfamily === subfamily;
            });
            
            // Create indexes
            await this.createIndexes(db);
            
        } finally {
            db.close();
        }
    }

    async createTableFromTSV(db) {
        // Read headers to create table
        const fileStream = createReadStream(this.inputPath);
        const rl = createInterface({ input: fileStream, crlfDelay: Infinity });
        
        const line = await rl[Symbol.asyncIterator]().next();
        const headers = line.value.split('\t');
        rl.close();
        
        const createSQL = `
            CREATE TABLE records (
                ${headers.map(header => `"${header.replace(/"/g, '""')}" TEXT`).join(',\n')}
            )
        `;
        
        await this.execAsync(db, createSQL);
    }

    async insertTSVData(db, filterFunction) {
        const fileStream = createReadStream(this.inputPath);
        const rl = createInterface({ input: fileStream, crlfDelay: Infinity });
        
        let headers = [];
        let lineCount = 0;
        let insertedCount = 0;
        
        // Prepare insert statement
        let stmt;
        
        await this.execAsync(db, 'BEGIN TRANSACTION');
        
        for await (const line of rl) {
            if (lineCount === 0) {
                headers = line.split('\t');
                const placeholders = headers.map(() => '?').join(',');
                stmt = db.prepare(`INSERT INTO records VALUES (${placeholders})`);
                lineCount++;
                continue;
            }
            
            const columns = line.split('\t');
            
            if (filterFunction(columns, headers)) {
                stmt.run(columns);
                insertedCount++;
            }
            
            lineCount++;
        }
        
        stmt.finalize();
        await this.execAsync(db, 'COMMIT');
        
        return insertedCount;
    }

    async analyzeSharedBINs(family, outputDir) {
        console.log(`    🔍 Analyzing shared BINs for ${family}...`);
        
        let sharedBins = [];
        
        if (this.isDatabase) {
            const db = new sqlite3.Database(this.inputPath, sqlite3.OPEN_READONLY);
            try {
                sharedBins = await this.queryAsync(db, `
                    SELECT 
                        bin_uri,
                        GROUP_CONCAT(DISTINCT COALESCE(subfamily, '${CONFIG.UNKNOWN_SUBFAMILY}')) as subfamilies,
                        COUNT(DISTINCT COALESCE(subfamily, '${CONFIG.UNKNOWN_SUBFAMILY}')) as subfamily_count
                    FROM records 
                    WHERE COALESCE(family, '${CONFIG.UNKNOWN_FAMILY}') = ? 
                        AND bin_uri IS NOT NULL 
                        AND bin_uri != 'None' 
                        AND bin_uri != ''
                    GROUP BY bin_uri 
                    HAVING subfamily_count > 1
                `, [family]);
            } finally {
                db.close();
            }
        } else {
            // Analyze from TSV (more complex, would need to scan file again)
            // For now, skip this for TSV input
            console.log(`    ⏭️  Skipping BIN analysis for TSV input`);
            return;
        }
        
        if (sharedBins.length > 0) {
            const csvPath = path.join(outputDir, 'subfamily_shared_BINs.csv');
            const csvContent = sharedBins.map(row => 
                `${row.bin_uri};${row.subfamilies};${family}`
            ).join('\n');
            
            fs.writeFileSync(csvPath, csvContent);
            console.log(`    📄 Saved ${sharedBins.length} shared BINs to ${csvPath}`);
        }
    }

    async createIndexes(db) {
        const indexes = [
            'CREATE INDEX IF NOT EXISTS idx_family ON records(family)',
            'CREATE INDEX IF NOT EXISTS idx_genus ON records(genus)',
            'CREATE INDEX IF NOT EXISTS idx_species ON records(species)',
            'CREATE INDEX IF NOT EXISTS idx_bin ON records(bin_uri)',
            'CREATE INDEX IF NOT EXISTS idx_processid ON records(processid)'
        ];
        
        for (const indexSQL of indexes) {
            try {
                await this.execAsync(db, indexSQL);
            } catch (err) {
                // Ignore errors for missing columns
            }
        }
    }

    sanitizeFileName(name) {
        if (!name || name === 'None' || name === '') {
            return 'Unknown';
        }
        return name.replace(/[<>:"/\\|?*]/g, '_').replace(/\s+/g, '_');
    }

    generateReport() {
        const reportPath = path.join(this.outputBase, 'splitting_report.txt');
        const report = [];
        
        report.push('BOLD Database Splitting Report');
        report.push('================================');
        report.push(`Input: ${this.inputPath}`);
        report.push(`Output: ${this.outputBase}`);
        report.push(`Threshold: ${CONFIG.FAMILY_SIZE_THRESHOLD} records`);
        report.push(`Generated: ${new Date().toISOString()}`);
        report.push('');
        
        report.push('Family Statistics:');
        report.push('-'.repeat(50));
        
        for (const [family, stats] of this.familyStats) {
            const strategy = stats.subfamilies.length > 0 ? 'subfamily' : 'family';
            report.push(`${family}: ${stats.count} records (split by ${strategy})`);
            
            if (stats.subfamilies.length > 0) {
                for (const subfam of stats.subfamilies) {
                    report.push(`  └─ ${subfam.subfamily}: ${subfam.count} records`);
                }
            }
        }
        
        report.push('');
        report.push('Output Structure:');
        report.push('-'.repeat(50));
        
        const totalFamilies = this.familyStats.size;
        const splitFamilies = Array.from(this.familyStats.values()).filter(s => s.subfamilies.length > 0).length;
        
        report.push(`Total families processed: ${totalFamilies}`);
        report.push(`Families split by subfamily: ${splitFamilies}`);
        report.push(`Regular family databases: ${totalFamilies - splitFamilies}`);
        
        fs.writeFileSync(reportPath, report.join('\n'));
        console.log(`\n📋 Report saved to: ${reportPath}`);
    }

    // Helper methods
    queryAsync(db, sql, params = []) {
        return new Promise((resolve, reject) => {
            db.all(sql, params, (err, rows) => {
                if (err) reject(err);
                else resolve(rows);
            });
        });
    }

    execAsync(db, sql, params = []) {
        return new Promise((resolve, reject) => {
            db.run(sql, params, function(err) {
                if (err) reject(err);
                else resolve(this);
            });
        });
    }
}

// Main execution
async function main() {
    if (process.argv.length < 3) {
        console.log(`
BOLD Family Database Splitter

Usage:
  node bold_splitter.js <input> [output_dir] [threshold]

Arguments:
  input        Path to input database (.db) or TSV file (.tsv)
  output_dir   Output directory (default: taxonomic_output)
  threshold    Family size threshold for subfamily splitting (default: 10000)

Examples:
  node bold_splitter.js bold.db
  node bold_splitter.js result_output.tsv my_output 5000
  node bold_splitter.js "C:\\GitHub\\bold-library-curation\\results\\_bags_sub-3\\bold.db"

Output Structure:
  output_dir/
  ├── Phylum1/
  │   ├── Class1/
  │   │   ├── Order1/
  │   │   │   ├── SmallFamily.db
  │   │   │   └── LargeFamily/
  │   │   │       ├── Subfamily1.db
  │   │   │       ├── Subfamily2.db
  │   │   │       └── subfamily_shared_BINs.csv
        `);
        process.exit(1);
    }

    const inputPath = process.argv[2];
    const outputDir = process.argv[3] || 'taxonomic_output';
    const threshold = parseInt(process.argv[4]) || 10000;
    
    // Validate input
    if (!fs.existsSync(inputPath)) {
        console.error(`❌ Input file not found: ${inputPath}`);
        process.exit(1);
    }
    
    // Set configuration
    CONFIG.FAMILY_SIZE_THRESHOLD = threshold;
    
    try {
        const splitter = new BoldFamilySplitter(inputPath, outputDir);
        await splitter.run();
    } catch (error) {
        console.error('❌ Processing failed:', error.message);
        process.exit(1);
    }
}

if (require.main === module) {
    main();
}

module.exports = { BoldFamilySplitter };