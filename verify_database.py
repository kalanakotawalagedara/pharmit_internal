"""Quick script to verify database structure"""
import sqlite3
from pathlib import Path

db_path = Path('Database/pharmacophores.db')

if not db_path.exists():
    print(f"ERROR: Database not found at {db_path}")
    exit(1)

conn = sqlite3.connect(db_path)
cursor = conn.cursor()

print("=" * 70)
print("DATABASE STRUCTURE")
print("=" * 70)

# List tables
cursor.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")
tables = [row[0] for row in cursor.fetchall()]
print(f"\nTables ({len(tables)}):")
for table in tables:
    print(f"  - {table}")

# Count rows in each table
print("\nRow counts:")
for table in tables:
    cursor.execute(f"SELECT COUNT(*) FROM {table}")
    count = cursor.fetchone()[0]
    print(f"  {table}: {count}")

# Show sample data from molecules table
print("\n" + "=" * 70)
print("MOLECULES TABLE SAMPLE")
print("=" * 70)
cursor.execute("SELECT * FROM molecules LIMIT 5")
for row in cursor.fetchall():
    print(row)

# Show sample triplets
print("\n" + "=" * 70)
print("TRIPLETS TABLE SAMPLE (first 5)")
print("=" * 70)
cursor.execute("SELECT * FROM triplets LIMIT 5")
for row in cursor.fetchall():
    print(row)

# Show sample features
print("\n" + "=" * 70)
print("FEATURES TABLE SAMPLE (first 10)")
print("=" * 70)
cursor.execute("SELECT * FROM features LIMIT 10")
for row in cursor.fetchall():
    print(row)

# Show triplet type distribution
print("\n" + "=" * 70)
print("TRIPLET TYPE DISTRIBUTION")
print("=" * 70)

# Get type names first
cursor.execute("SELECT type_id, type_name FROM pharma_types ORDER BY type_id")
type_names = {row[0]: row[1] for row in cursor.fetchall()}

cursor.execute("""
    SELECT type1_id, type2_id, type3_id, COUNT(*) as count 
    FROM triplets 
    GROUP BY type1_id, type2_id, type3_id 
    ORDER BY count DESC 
    LIMIT 10
""")

for row in cursor.fetchall():
    t1 = type_names.get(row[0], f"Unknown({row[0]})")
    t2 = type_names.get(row[1], f"Unknown({row[1]})")
    t3 = type_names.get(row[2], f"Unknown({row[2]})")
    print(f"  ({t1}, {t2}, {t3}): {row[3]}")

# Show indices
print("\n" + "=" * 70)
print("INDICES")
print("=" * 70)
cursor.execute("SELECT name, tbl_name FROM sqlite_master WHERE type='index' ORDER BY name")
for row in cursor.fetchall():
    if row[0].startswith('sqlite_'):
        continue
    print(f"  {row[0]} on {row[1]}")

conn.close()
print("\n✓ Database verification complete!")
