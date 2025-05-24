const fs = require('fs');
const path = require('path');
const sqlite3 = require('sqlite3').verbose();
const { parseString } = require('xml2js');
const { promisify } = require('util');

const DATA_FOLDER = path.join(__dirname, 'xml_files');
const OUTPUT_FOLDER = path.join(__dirname, 'db_files');

const parseStringAsync = promisify(parseString);

// Recursively find all .xml files in a directory
function findXmlFiles(dir) {
    let xmlFiles = [];
    const entries = fs.readdirSync(dir, { withFileTypes: true });

    for (const entry of entries) {
        const fullPath = path.join(dir, entry.name);
        if (entry.isDirectory()) {
            xmlFiles = xmlFiles.concat(findXmlFiles(fullPath));
        } else if (entry.isFile() && entry.name.endsWith('.xml')) {
            xmlFiles.push(fullPath);
        }
    }

    return xmlFiles;
}

// Normalize field values
function normalizeValue(value, col) {
    if (value === undefined || value === null || value === 'None') return null;

    // Special formatting for coord field
    if (col === 'coord' && Array.isArray(value) && value.length === 2) {
        return `${value[0]}, ${value[1]}`;
    }

    if (Array.isArray(value)) {
        const flat = value.map(v => {
            if (typeof v === 'object' && v !== null && 'id' in v) return v.id;
            return typeof v === 'object' ? JSON.stringify(v) : v;
        });
        return [...new Set(flat)].filter(Boolean).join('; ') || null;
    }

    if (typeof value === 'object' && value !== null) {
        if ('id' in value) return value.id;
        return JSON.stringify(value);
    }

    if (typeof value === 'string' && value.trim().toLowerCase() === 'none') {
        return null;
    }

    return value;
}

// Override and prioritize some field mappings
const fieldOverrides = {
    recordid: (r) => r.recordid || r.record_id || null,
    country: (r) => r.country_ocean || r.country || null,
    COLLECTION_DATE: (r) => r.collection_date || r.collection_date_start || null
};

function getValue(record, col) {
    if (col === 'url') {
        return record.processid ? `https://portal.boldsystems.org/record/${record.processid}` : null;
    }

    const rawValue = fieldOverrides[col] ? fieldOverrides[col](record) : record[col];
    return normalizeValue(rawValue, col);
}

function initDatabaseAndLoadData(xmlPath, database) {
    return new Promise((resolve, reject) => {
        fs.readFile(xmlPath, 'utf8', async (readErr, xmlData) => {
            if (readErr) return reject(`Error reading XML: ${readErr}`);

            let result;
            try {
                result = await parseStringAsync(xmlData, { explicitArray: false });
            } catch (parseErr) {
                return reject(`Error parsing XML: ${parseErr}`);
            }

            if (!result.records || !result.records.record) {
                console.warn(`⚠️ No records found in ${xmlPath}. Skipping.`);
                return resolve(); // continue with next file
            }

            const records = Array.isArray(result.records.record)
                ? result.records.record
                : [result.records.record];

            const columns = [
                'url','keep','ranking','BAGS','status','recordid','taxonid','processid','sampleid','fieldid',
                'museumid','record_id','specimenid','processid_minted_date','bin_uri','bin_created_date',
                'collection_code','inst','taxid','kingdom','phylum','class','order','family','subfamily','tribe',
                'genus','species','subspecies','species_reference','identification','identification_method',
                'identification_rank','identified_by','identifier_email','taxonomy_notes','sex','reproduction',
                'life_stage','short_note','notes','voucher_type','tissue_type','specimen_linkout',
                'associated_specimens','associated_taxa','collection_date','collection_date_start', 'collection_date_end',
                'collection_date_accuracy','collection_event_id','collection_time','collection_notes','geoid',
                'country_ocean','country_iso','province_state','region','sector','site','site_code','coord',
                'coord_accuracy','coord_source','elev','elev_accuracy','depth','depth_accuracy','habitat',
                'sampling_protocol','nuc','nuc_basecount','insdc_acs','funding_src','marker_code',
                'primers_forward','primers_reverse','sequence_run_site','sequence_upload_date','recordset_code_arr',
                'extrainfo','country','collection_note','associated_specimen','gb_acs','nucraw','SPECIES_ID',
                'TYPE_SPECIMEN','SEQ_QUALITY','HAS_IMAGE','HAS_COLLECTOR','IDENTIFIER','ID_METHOD','INSTITUTION',
                'PUBLIC_VOUCHER','MUSEUM_ID','curator_notes','collectors'
            ];

            const createTableSQL = `CREATE TABLE IF NOT EXISTS records (
                ${columns.map(col => `"${col}" TEXT`).join(',\n')}
            )`;

            database.serialize(() => {
                database.run("DROP TABLE IF EXISTS records");
                database.run(createTableSQL);

                const placeholders = columns.map(() => '?').join(',');
                const insertStmt = database.prepare(`INSERT OR REPLACE INTO records (${columns.map(c => `"${c}"`).join(',')}) VALUES (${placeholders})`);

                for (const record of records) {
                    const values = columns.map(col => getValue(record, col));
                    insertStmt.run(values);
                }

                insertStmt.finalize();
                resolve();
            });
        });
    });
}

// Ensure output directory exists
fs.mkdirSync(OUTPUT_FOLDER, { recursive: true });

// Process XML files
const xmlFiles = findXmlFiles(DATA_FOLDER);

xmlFiles.forEach(xmlPath => {
    const relativePath = path.relative(DATA_FOLDER, xmlPath);
    const dbRelative = relativePath.replace(/\.xml$/, '.db');
    const dbPath = path.join(OUTPUT_FOLDER, dbRelative);
    const dbDir = path.dirname(dbPath);

    fs.mkdirSync(dbDir, { recursive: true });

    const db = new sqlite3.Database(dbPath, err => {
        if (err) return console.error(`DB error for ${dbPath}: ${err.message}`);

        initDatabaseAndLoadData(xmlPath, db)
            .then(() => {
                console.log(`✅ Created: ${dbPath}`);
                db.close();
            })
            .catch(err => {
                console.error(`❌ Failed for ${xmlPath}:`, err);
                db.close();
            });
    });
});