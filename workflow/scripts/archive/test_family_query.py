#!/usr/bin/env python3
"""
Test script to verify the SQL syntax fix
"""

import sqlite3
import sys

def test_family_query(db_path):
    """Test the family counting query"""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Test the corrected query
        query = """
        SELECT 
            COALESCE(kingdom, 'Unknown') as kingdom,
            COALESCE(phylum, 'Unknown') as phylum,
            COALESCE("class", 'Unknown') as class,
            COALESCE("order", 'Unknown') as "order",
            COALESCE(family, 'Unknown') as family,
            COUNT(*) as record_count,
            COUNT(DISTINCT subfamily) as subfamily_count
        FROM (
            SELECT DISTINCT recordid, kingdom, phylum, "class", "order", family, subfamily
            FROM bold
            WHERE family IS NOT NULL AND family != ''
        )
        GROUP BY kingdom, phylum, "class", "order", family
        HAVING record_count > 0
        ORDER BY record_count DESC
        LIMIT 5
        """
        
        cursor.execute(query)
        families = cursor.fetchall()
        
        print(f"✓ Query executed successfully!")
        print(f"Found {len(families)} families (showing top 5):")
        
        for i, row in enumerate(families):
            kingdom, phylum, class_name, order, family, count, subfamily_count = row
            print(f"  {i+1}. {kingdom}/{phylum}/{class_name}/{order}/{family} - {count} records, {subfamily_count} subfamilies")
        
        conn.close()
        return True
        
    except Exception as e:
        print(f"✗ Query failed: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python test_family_query.py <database_path>")
        sys.exit(1)
    
    db_path = sys.argv[1]
    print(f"Testing family query on: {db_path}")
    
    if test_family_query(db_path):
        print("\n✓ SQL syntax fix verified - prepare_family_batches.py should work now")
    else:
        print("\n✗ There may be other issues with the database or query")
