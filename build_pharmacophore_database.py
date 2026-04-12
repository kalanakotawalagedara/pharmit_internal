#!/usr/bin/env python3
"""
Pharmacophore Database Builder
Builds searchable SQLite database with triangle decomposition and spatial indexing
Based on Pharmit's database architecture for recursive backtracking search

Usage:
    python build_pharmacophore_database.py --json-dir Output --output Database
"""

import sqlite3
import json
import itertools
import logging
import sys
from pathlib import Path
from typing import Dict, Optional
from datetime import datetime
import numpy as np

# Logger will be configured in main() based on command-line arguments
logger = logging.getLogger(__name__)


# Pharmacophore type definitions (matching Pharmit's order)
PHARMACOPHORE_TYPES = {
    'Aromatic': 0,
    'HydrogenDonor': 1,
    'HydrogenAcceptor': 2,
    'PositiveIon': 3,
    'NegativeIon': 4,
    'Hydrophobic': 5
}


class PharmacophoreDatabase:
    """
    Build searchable pharmacophore database with triangle decomposition
    Implements SQL-based approximation of Pharmit's KDB-tree indexing
    """
    
    def __init__(self, db_path: Path):
        self.db_path = db_path
        self.conn = None
        self.stats = {
            'molecules_processed': 0,
            'molecules_failed': 0,
            'molecules_skipped': 0,
            'features_stored': 0,
            'triplets_generated': 0,
            'vectors_stored': 0,
            'errors': []
        }
    
    def connect(self):
        """Initialize database connection"""
        self.conn = sqlite3.connect(str(self.db_path))
        self.conn.row_factory = sqlite3.Row
        logger.info(f"Connected to database: {self.db_path}")
    
    def create_schema(self):
        """Create database schema with indices for fast search"""
        if not self.conn:
            raise RuntimeError("Database not connected")
        
        cursor = self.conn.cursor()
        
        logger.info("Creating database schema...")
        
        # Pharmacophore type definitions
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS pharma_types (
                type_id INTEGER PRIMARY KEY,
                type_name TEXT UNIQUE NOT NULL
            )
        ''')
        
        # Insert pharmacophore types
        for type_name, type_id in PHARMACOPHORE_TYPES.items():
            cursor.execute(
                'INSERT OR IGNORE INTO pharma_types (type_id, type_name) VALUES (?, ?)',
                (type_id, type_name)
            )
        
        # Molecule registry (PDB → mol_id mapping)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS molecules (
                mol_id INTEGER PRIMARY KEY AUTOINCREMENT,
                pdb_id TEXT NOT NULL,
                ligand_name TEXT NOT NULL,
                chain TEXT,
                num_features INTEGER NOT NULL,
                json_path TEXT NOT NULL,
                created_date TIMESTAMP DEFAULT CURRENT_TIMESTAMP,
                UNIQUE(pdb_id, ligand_name)
            )
        ''')
        
        # Triangle decomposition (triplets)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS triplets (
                triplet_id INTEGER PRIMARY KEY AUTOINCREMENT,
                mol_id INTEGER NOT NULL,
                type1_id INTEGER NOT NULL,
                type2_id INTEGER NOT NULL,
                type3_id INTEGER NOT NULL,
                dist12 REAL NOT NULL,
                dist23 REAL NOT NULL,
                dist31 REAL NOT NULL,
                centroid_x REAL,
                centroid_y REAL,
                centroid_z REAL,
                point1_idx INTEGER NOT NULL,
                point2_idx INTEGER NOT NULL,
                point3_idx INTEGER NOT NULL,
                FOREIGN KEY(mol_id) REFERENCES molecules(mol_id),
                FOREIGN KEY(type1_id) REFERENCES pharma_types(type_id),
                FOREIGN KEY(type2_id) REFERENCES pharma_types(type_id),
                FOREIGN KEY(type3_id) REFERENCES pharma_types(type_id)
            )
        ''')
        
        # H-bond vectors (directional information)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS vectors (
                vector_id INTEGER PRIMARY KEY AUTOINCREMENT,
                triplet_id INTEGER NOT NULL,
                point_idx INTEGER NOT NULL,
                vx REAL NOT NULL,
                vy REAL NOT NULL,
                vz REAL NOT NULL,
                FOREIGN KEY(triplet_id) REFERENCES triplets(triplet_id)
            )
        ''')
        
        # Pharmacophore features (all points with coordinates for RMSD calculation)
        cursor.execute('''
            CREATE TABLE IF NOT EXISTS features (
                feature_id INTEGER PRIMARY KEY AUTOINCREMENT,
                mol_id INTEGER NOT NULL,
                point_idx INTEGER NOT NULL,
                type_id INTEGER NOT NULL,
                x REAL NOT NULL,
                y REAL NOT NULL,
                z REAL NOT NULL,
                vector_x REAL,
                vector_y REAL,
                vector_z REAL,
                FOREIGN KEY(mol_id) REFERENCES molecules(mol_id),
                FOREIGN KEY(type_id) REFERENCES pharma_types(type_id)
            )
        ''')
        
        logger.info("Creating indices for fast search...")
        
        # Critical indices for recursive backtracking search
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_triplet_types ON triplets(type1_id, type2_id, type3_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_dist12 ON triplets(dist12)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_dist23 ON triplets(dist23)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_dist31 ON triplets(dist31)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_mol_lookup ON molecules(pdb_id, ligand_name)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_spatial ON triplets(centroid_x, centroid_y, centroid_z)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_mol_id ON triplets(mol_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_vector_triplet ON vectors(triplet_id)')
        cursor.execute('CREATE INDEX IF NOT EXISTS idx_mol_features ON features(mol_id)')
        
        self.conn.commit()
        logger.info("Schema created successfully")
    
    def calculate_distance(self, p1: Dict, p2: Dict) -> float:
        """Calculate Euclidean distance between two points"""
        dx = p1['x'] - p2['x']
        dy = p1['y'] - p2['y']
        dz = p1['z'] - p2['z']
        return float(np.sqrt(dx*dx + dy*dy + dz*dz))
    
    def add_molecule(self, json_path: Path) -> bool:
        """
        Add molecule from JSON and decompose into triplets
        Returns: True if successful, False if failed
        """
        if not self.conn:
            raise RuntimeError("Database not connected")
        
        try:
            # Load JSON
            with open(json_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            pdb_id = data.get('pdb_id', '').upper()
            ligand_name = data.get('ligand_name', '')
            chain = data.get('ligand_chain', '')
            features = data.get('features', [])
            
            # Validate required fields
            if not pdb_id or not ligand_name:
                logger.error(f"Missing PDB ID or ligand name in {json_path}")
                self.stats['molecules_failed'] += 1
                self.stats['errors'].append({
                    'file': str(json_path),
                    'error': 'Missing PDB ID or ligand name'
                })
                return False
            
            # Need at least 3 features for triplet generation
            if len(features) < 3:
                logger.warning(f"{pdb_id}_{ligand_name}: Only {len(features)} features, skipping")
                self.stats['molecules_skipped'] += 1
                return True
            
            logger.info(f"Processing {pdb_id}_{ligand_name} ({len(features)} features)")
            
            cursor = self.conn.cursor()
            
            # Insert molecule (skip if duplicate)
            try:
                cursor.execute('''
                    INSERT INTO molecules (pdb_id, ligand_name, chain, num_features, json_path)
                    VALUES (?, ?, ?, ?, ?)
                ''', (pdb_id, ligand_name, chain, len(features), str(json_path)))
                mol_id = cursor.lastrowid
            except sqlite3.IntegrityError:
                logger.warning(f"{pdb_id}_{ligand_name}: Already exists, skipping")
                self.stats['molecules_skipped'] += 1
                return True
            
            # Store all pharmacophore features with coordinates
            feature_count = 0
            for point_idx, feature in enumerate(features):
                # Get type ID
                if feature['type'] not in PHARMACOPHORE_TYPES:
                    logger.warning(f"Unknown pharmacophore type '{feature['type']}', skipping feature")
                    continue
                
                type_id = PHARMACOPHORE_TYPES[feature['type']]
                
                # Extract coordinates
                x = feature.get('x', 0.0)
                y = feature.get('y', 0.0)
                z = feature.get('z', 0.0)
                
                # Extract vector if present
                vector_x, vector_y, vector_z = None, None, None
                if 'vector' in feature and feature['vector']:
                    vec = feature['vector']
                    vector_x = vec.get('x')
                    vector_y = vec.get('y')
                    vector_z = vec.get('z')
                
                # Insert feature
                cursor.execute('''
                    INSERT INTO features (
                        mol_id, point_idx, type_id, x, y, z,
                        vector_x, vector_y, vector_z
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (mol_id, point_idx, type_id, x, y, z, vector_x, vector_y, vector_z))
                
                feature_count += 1
            
            # Generate all C(N,3) triplets
            triplet_count = 0
            vector_count = 0
            
            for i, j, k in itertools.combinations(range(len(features)), 3):
                p1, p2, p3 = features[i], features[j], features[k]
                
                # Validate feature types
                if p1['type'] not in PHARMACOPHORE_TYPES or \
                   p2['type'] not in PHARMACOPHORE_TYPES or \
                   p3['type'] not in PHARMACOPHORE_TYPES:
                    logger.warning(f"Unknown pharmacophore type in triplet, skipping")
                    continue
                
                # Get type IDs and sort for canonical representation
                type_ids = sorted([
                    PHARMACOPHORE_TYPES[p1['type']],
                    PHARMACOPHORE_TYPES[p2['type']],
                    PHARMACOPHORE_TYPES[p3['type']]
                ])
                
                # Calculate triangle edge distances
                dist12 = self.calculate_distance(p1, p2)
                dist23 = self.calculate_distance(p2, p3)
                dist31 = self.calculate_distance(p3, p1)
                
                # Calculate centroid for spatial indexing
                centroid_x = (p1['x'] + p2['x'] + p3['x']) / 3.0
                centroid_y = (p1['y'] + p2['y'] + p3['y']) / 3.0
                centroid_z = (p1['z'] + p2['z'] + p3['z']) / 3.0
                
                # Insert triplet
                cursor.execute('''
                    INSERT INTO triplets (
                        mol_id, type1_id, type2_id, type3_id,
                        dist12, dist23, dist31,
                        centroid_x, centroid_y, centroid_z,
                        point1_idx, point2_idx, point3_idx
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ''', (
                    mol_id, type_ids[0], type_ids[1], type_ids[2],
                    dist12, dist23, dist31,
                    centroid_x, centroid_y, centroid_z,
                    i, j, k
                ))
                
                triplet_id = cursor.lastrowid
                triplet_count += 1
                
                # Store H-bond vectors
                for idx, point in enumerate([p1, p2, p3]):
                    if 'vector' in point and point['vector']:
                        vec = point['vector']
                        cursor.execute('''
                            INSERT INTO vectors (triplet_id, point_idx, vx, vy, vz)
                            VALUES (?, ?, ?, ?, ?)
                        ''', (triplet_id, idx, vec['x'], vec['y'], vec['z']))
                        vector_count += 1
            
            self.conn.commit()
            
            logger.info(f"  Stored {feature_count} features, generated {triplet_count} triplets, {vector_count} vectors")
            self.stats['molecules_processed'] += 1
            self.stats['features_stored'] += feature_count
            self.stats['triplets_generated'] += triplet_count
            self.stats['vectors_stored'] += vector_count
            
            return True
            
        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON in {json_path}: {e}")
            self.stats['molecules_failed'] += 1
            self.stats['errors'].append({
                'file': str(json_path),
                'error': f'JSON decode error: {str(e)}'
            })
            return False
        except KeyError as e:
            logger.error(f"Missing required field in {json_path}: {e}")
            self.stats['molecules_failed'] += 1
            self.stats['errors'].append({
                'file': str(json_path),
                'error': f'Missing field: {str(e)}'
            })
            return False
        except Exception as e:
            logger.error(f"Unexpected error processing {json_path}: {e}")
            self.stats['molecules_failed'] += 1
            self.stats['errors'].append({
                'file': str(json_path),
                'error': str(e)
            })
            return False
    
    def build_from_directory(self, json_dir: Path):
        """Build database from directory of JSON files"""
        json_files = sorted(json_dir.glob('*.json'))
        
        if not json_files:
            logger.error(f"No JSON files found in {json_dir}")
            return
        
        logger.info(f"Found {len(json_files)} JSON files")
        logger.info("="*70)
        
        for i, json_file in enumerate(json_files, 1):
            if i % 100 == 0:
                logger.info(f"Progress: {i}/{len(json_files)} files processed")
            
            self.add_molecule(json_file)
        
        logger.info("="*70)
        logger.info("Database build complete!")
    
    def export_molecule_index(self, output_path: Path):
        """Export molecule index as CSV for quick PDB lookup"""
        if not self.conn:
            raise RuntimeError("Database not connected")
        
        cursor = self.conn.cursor()
        cursor.execute('''
            SELECT mol_id, pdb_id, ligand_name, chain, num_features, json_path
            FROM molecules
            ORDER BY mol_id
        ''')
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write("mol_id,pdb_id,ligand_name,chain,num_features,json_path\n")
            for row in cursor.fetchall():
                f.write(f"{row['mol_id']},{row['pdb_id']},{row['ligand_name']},"
                       f"{row['chain']},{row['num_features']},{row['json_path']}\n")
        
        logger.info(f"Molecule index exported to {output_path}")
    
    def print_statistics(self):
        """Print database statistics"""
        if not self.conn:
            raise RuntimeError("Database not connected")
        
        cursor = self.conn.cursor()
        
        # Get counts
        cursor.execute('SELECT COUNT(*) as count FROM molecules')
        num_mols = cursor.fetchone()['count']
        
        cursor.execute('SELECT COUNT(*) as count FROM features')
        num_features = cursor.fetchone()['count']
        
        cursor.execute('SELECT COUNT(*) as count FROM triplets')
        num_triplets = cursor.fetchone()['count']
        
        cursor.execute('SELECT COUNT(*) as count FROM vectors')
        num_vectors = cursor.fetchone()['count']
        
        cursor.execute('SELECT AVG(num_features) as avg FROM molecules')
        avg_features = cursor.fetchone()['avg'] or 0
        
        # Get triplet type distribution
        cursor.execute('''
            SELECT type1_id, type2_id, type3_id, COUNT(*) as count
            FROM triplets
            GROUP BY type1_id, type2_id, type3_id
            ORDER BY count DESC
            LIMIT 10
        ''')
        top_triplet_types = cursor.fetchall()
        
        print("\n" + "="*70)
        print("DATABASE STATISTICS")
        print("="*70)
        print(f"Molecules processed:    {self.stats['molecules_processed']:,}")
        print(f"Molecules failed:       {self.stats['molecules_failed']:,}")
        print(f"Molecules skipped:      {self.stats['molecules_skipped']:,}")
        print(f"Total in database:      {num_mols:,}")
        print(f"-" * 70)
        print(f"Total features:         {num_features:,}")
        print(f"Total triplets:         {num_triplets:,}")
        print(f"Total vectors:          {num_vectors:,}")
        print(f"Avg features/molecule:  {avg_features:.1f}")
        print(f"Avg triplets/molecule:  {num_triplets/num_mols if num_mols > 0 else 0:.1f}")
        print(f"-" * 70)
        
        if top_triplet_types:
            print("Top 10 triplet type combinations:")
            type_names = {v: k for k, v in PHARMACOPHORE_TYPES.items()}
            for row in top_triplet_types:
                t1 = type_names.get(row['type1_id'], '?')
                t2 = type_names.get(row['type2_id'], '?')
                t3 = type_names.get(row['type3_id'], '?')
                print(f"  ({t1}, {t2}, {t3}): {row['count']:,}")
        
        print("="*70 + "\n")
    
    def save_build_log(self, log_path: Path):
        """Save build log with statistics and errors"""
        log_data = {
            'build_date': datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3] + 'Z',
            'database_path': str(self.db_path),
            'statistics': self.stats,
            'pharmacophore_types': PHARMACOPHORE_TYPES
        }
        
        with open(log_path, 'w', encoding='utf-8') as f:
            json.dump(log_data, f, indent=2)
        
        logger.info(f"Build log saved to {log_path}")
    
    def close(self):
        """Close database connection"""
        if self.conn:
            self.conn.close()
            logger.info("Database connection closed")


def main():
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Build pharmacophore database with triangle decomposition',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Build database from JSON directory
  python build_pharmacophore_database.py --json-dir Output --output Database

  # Specify custom output paths
  python build_pharmacophore_database.py \\
    --json-dir Output \\
    --output Database \\
    --database-name pharmacophores.db \\
    --index-name molecule_index.csv
        '''
    )
    
    parser.add_argument('--json-dir', type=str, required=True,
                       help='Directory containing pharmacophore JSON files')
    parser.add_argument('--output', type=str, default='Database',
                       help='Output directory for database files')
    parser.add_argument('--database-name', type=str, default='pharmacophores.db',
                       help='Database filename')
    parser.add_argument('--index-name', type=str, default='molecule_index.csv',
                       help='Molecule index filename')
    parser.add_argument('--log-name', type=str, default='build_log.json',
                       help='Build log filename')
    parser.add_argument('--log-level', type=str, default='INFO',
                       choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                       help='Logging level')
    
    args = parser.parse_args()
    
    # Setup logging
    numeric_level = getattr(logging, args.log_level.upper(), None)
    logging.basicConfig(
        level=numeric_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler()
        ]
    )
    
    # Setup paths
    json_dir = Path(args.json_dir)
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    db_path = output_dir / args.database_name
    index_path = output_dir / args.index_name
    log_path = output_dir / args.log_name
    
    if not json_dir.exists():
        logger.error(f"JSON directory not found: {json_dir}")
        return 1
    
    # Build database
    logger.info("="*70)
    logger.info("PHARMACOPHORE DATABASE BUILDER")
    logger.info("="*70)
    logger.info(f"Input directory:  {json_dir}")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Database:         {db_path}")
    logger.info("="*70 + "\n")
    
    db = PharmacophoreDatabase(db_path)
    
    try:
        db.connect()
        db.create_schema()
        db.build_from_directory(json_dir)
        db.print_statistics()
        db.export_molecule_index(index_path)
        db.save_build_log(log_path)
        
    except Exception as e:
        logger.error(f"Fatal error: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    finally:
        db.close()
    
    logger.info("✓ Database build complete!")
    logger.info(f"  Database: {db_path}")
    logger.info(f"  Index:    {index_path}")
    logger.info(f"  Log:      {log_path}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
