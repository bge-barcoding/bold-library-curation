# Database Wrangling Scripts

This directory contains two Node.js scripts designed to process and organize BOLD (Barcode of Life Data Systems) biological data from XML format into SQLite databases with taxonomic organization.

## Overview

The database wrangling workflow consists of two main steps:
1. **XML to Database Conversion** (`build_dbs.js`) - Converts XML files containing biological records into SQLite databases
2. **Taxonomic Organization** (`sort_existing_dbs.js`) - Sorts and organizes the generated databases by taxonomic hierarchy

## Scripts

### 1. build_dbs.js

Converts XML files containing biological specimen records into SQLite databases.

#### Purpose
- Recursively scans for XML files in the `xml_files` directory
- Parses each XML file and extracts biological record data
- Creates corresponding SQLite databases in the `db_files` directory
- Maintains directory structure from input to output

#### Key Features
- **Recursive XML Discovery**: Automatically finds all `.xml` files in subdirectories
- **Data Normalization**: Handles various data types and formats consistently
- **Field Mapping**: Maps XML elements to standardized database columns
- **URL Generation**: Creates BOLD portal URLs using process IDs
- **Error Handling**: Continues processing even if individual files fail

#### Database Schema
Creates a `records` table with 80+ columns including:
- **Identifiers**: `recordid`, `processid`, `sampleid`, `fieldid`, `museumid`
- **Taxonomy**: `kingdom`, `phylum`, `class`, `order`, `family`, `genus`, `species`
- **Collection Data**: `collection_date`, `country`, `province_state`, `coord`, `habitat`
- **Sequence Data**: `nuc`, `marker_code`, `primers_forward`, `primers_reverse`
- **Metadata**: `status`, `ranking`, `notes`, `voucher_type`, `tissue_type`

#### Field Processing
- **Coordinate Formatting**: Converts coordinate arrays to "lat, lon" format
- **Array Flattening**: Joins array values with semicolons
- **Object Handling**: Extracts IDs from nested objects or converts to JSON
- **Null Standardization**: Converts "None" and empty values to NULL

### 2. sort_existing_dbs.js

Organizes SQLite databases into a hierarchical directory structure based on taxonomic classification.

#### Purpose
- Reads existing SQLite databases from `db_files` directory
- Extracts taxonomic information from database records
- Reorganizes databases into taxonomic hierarchy in `sorted_dbs` directory
- Creates family-level folders when appropriate

#### Organization Logic
The script creates directory structures following taxonomic hierarchy:
```
sorted_dbs/
├── [Phylum]/
│   ├── [Class]/
│   │   ├── [Order]/
│   │   │   ├── [Family]/  (if family folder exists in input)
│   │   │   └── database.db
│   │   └── another_database.db
```

#### Key Features
- **Taxonomy Extraction**: Reads taxonomic data from the first record in each database
- **Intelligent Folder Detection**: Checks if family-level organization exists in input
- **Directory Creation**: Automatically creates nested taxonomic directories
- **File Safety**: Uses copy operations to preserve original files
- **Missing Data Handling**: Uses "Unknown_[level]" for missing taxonomic information

## Usage

### Prerequisites
```bash
npm install sqlite3 xml2js
```

### Directory Structure
Before running the scripts, ensure this directory structure:
```
database_wrangling/
├── xml_files/          # Input XML files (created by you)
├── db_files/           # Generated SQLite databases (created by build_dbs.js)
├── sorted_dbs/         # Taxonomically organized databases (created by sort_existing_dbs.js)
├── build_dbs.js
├── sort_existing_dbs.js
└── README.md
```

### Step 1: Convert XML to Databases
```bash
node build_dbs.js
```
- Place your XML files in the `xml_files` directory
- Script will create corresponding `.db` files in `db_files`
- Maintains subdirectory structure from input

### Step 2: Organize by Taxonomy
```bash
node sort_existing_dbs.js
```
- Processes all `.db` files from `db_files`
- Creates taxonomically organized structure in `sorted_dbs`
- Files are copied, not moved (originals preserved)

## Data Flow

```
XML Files (xml_files/) 
    ↓ build_dbs.js
SQLite Databases (db_files/)
    ↓ sort_existing_dbs.js  
Organized Databases (sorted_dbs/Phylum/Class/Order/[Family]/)
```

## Error Handling

Both scripts include robust error handling:
- **build_dbs.js**: Continues processing if individual XML files fail to parse
- **sort_existing_dbs.js**: Skips databases that can't be read or lack taxonomic data
- Console output provides clear success/failure indicators

## Output Examples

### build_dbs.js Output
```
✅ Created: C:\path\to\db_files\arthropoda\insects\beetles\beetle_data.db
✅ Created: C:\path\to\db_files\chordata\mammals\carnivores\carnivore_data.db
⚠️ No records found in empty_file.xml. Skipping.
```

### sort_existing_dbs.js Output
```
📦 Moved beetle_data.db → Arthropoda/Insecta/Coleoptera/Scarabaeidae/beetle_data.db
📦 Moved carnivore_data.db → Chordata/Mammalia/Carnivora/carnivore_data.db
⚠️ Skipping corrupt_file.db: Missing taxonomy
```

## Technical Notes

- **Database Format**: SQLite 3
- **Character Encoding**: UTF-8
- **Array Separators**: Semicolons (`;`)
- **Coordinate Format**: "latitude, longitude"
- **URL Pattern**: `https://portal.boldsystems.org/record/{processid}`

## Dependencies

- `sqlite3`: SQLite database operations
- `xml2js`: XML parsing and conversion
- `fs`: File system operations (Node.js built-in)
- `path`: Path manipulation utilities (Node.js built-in)

## Troubleshooting

### Common Issues
1. **Missing xml_files directory**: Create the directory and add your XML files
2. **SQLite errors**: Ensure proper read/write permissions on directories
3. **XML parsing failures**: Check XML file formatting and structure
4. **Memory issues**: For very large XML files, consider processing in batches

### Performance Tips
- Process XML files in smaller batches for better memory management
- Ensure adequate disk space for database files (typically 2-5x XML file size)
- Use SSD storage for faster database operations

## License

This code is part of the bold-library-curation project. See the main project README for licensing information.