const fs = require('fs');
const path = require('path');
const sqlite3 = require('sqlite3').verbose();

const INPUT_FOLDER = path.join(__dirname, 'db_files');
const OUTPUT_FOLDER = path.join(__dirname, 'sorted_dbs');

// Create output root folder if missing
if (!fs.existsSync(OUTPUT_FOLDER)) {
    fs.mkdirSync(OUTPUT_FOLDER);
}

// Recursively find all .db files
function findAllDbFiles(dir, collected = []) {
    const entries = fs.readdirSync(dir, { withFileTypes: true });
    for (const entry of entries) {
        const fullPath = path.join(dir, entry.name);
        if (entry.isDirectory()) {
            findAllDbFiles(fullPath, collected);
        } else if (entry.isFile() && entry.name.endsWith('.db')) {
            collected.push(fullPath);
        }
    }
    return collected;
}

// Get immediate family folders in INPUT
function getImmediateFamilyFolders(inputRoot) {
    const familyFolders = new Set();
    const entries = fs.readdirSync(inputRoot, { withFileTypes: true });
    for (const entry of entries) {
        if (entry.isDirectory()) {
            const familyPath = path.join(inputRoot, entry.name);
            const inner = fs.readdirSync(familyPath, { withFileTypes: true });
            const hasDb = inner.some(e => e.isFile() && e.name.endsWith('.db'));
            if (hasDb) {
                familyFolders.add(entry.name);
            }
        }
    }
    return familyFolders;
}

// Read taxonomy from one record
function getTaxonomy(filePath) {
    return new Promise((resolve) => {
        const db = new sqlite3.Database(filePath, sqlite3.OPEN_READONLY, err => {
            if (err) return resolve(null);
        });

        db.get("SELECT phylum, class, `order`, family FROM records LIMIT 1", (err, row) => {
            db.close();
            if (err || !row) return resolve(null);
            resolve({
                phylum: row.phylum || 'Unknown_phylum',
                class: row.class || 'Unknown_class',
                order: row.order || 'Unknown_order',
                family: row.family || 'Unknown_family',
            });
        });
    });
}

// MAIN
(async () => {
    const allDbFiles = findAllDbFiles(INPUT_FOLDER);
    const knownFamilyFolders = getImmediateFamilyFolders(INPUT_FOLDER);

    for (const filePath of allDbFiles) {
        const fileName = path.basename(filePath);
        const tax = await getTaxonomy(filePath);

        if (!tax) {
            console.warn(`⚠️ Skipping ${fileName}: Missing taxonomy`);
            continue;
        }

        const { phylum, class: cls, order, family } = tax;

        let destDir;
        if (knownFamilyFolders.has(family)) {
            destDir = path.join(OUTPUT_FOLDER, phylum, cls, order, family);
        } else {
            destDir = path.join(OUTPUT_FOLDER, phylum, cls, order);
        }

        fs.mkdirSync(destDir, { recursive: true });
        const destPath = path.join(destDir, fileName);
        fs.copyFileSync(filePath, destPath);

        console.log(`📦 Moved ${fileName} → ${path.relative(OUTPUT_FOLDER, destPath)}`);
    }
})();
