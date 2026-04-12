#!/usr/bin/env python3
"""
Pharmacophore Database Search Engine
Production-ready implementation with recursive backtracking and Kabsch RMSD calculation
"""

import sys
import json
import logging
import sqlite3
import argparse
import itertools
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Set
from collections import defaultdict
from datetime import datetime

try:
    import numpy as np
    import pandas as pd
except ImportError as e:
    print(f"ERROR: Missing required package: {e}")
    print("Install with: pip install numpy pandas")
    sys.exit(1)

# Configure logging
logger = logging.getLogger(__name__)

# Pharmacophore type definitions (matching database schema)
PHARMACOPHORE_TYPES = {
    'Aromatic': 0,
    'HydrogenDonor': 1,
    'HydrogenAcceptor': 2,
    'PositiveIon': 3,
    'NegativeIon': 4,
    'Hydrophobic': 5
}

TYPE_ID_TO_NAME = {v: k for k, v in PHARMACOPHORE_TYPES.items()}


class DatabaseConnectionError(Exception):
    """Raised when database connection fails"""
    pass


class QueryParseError(Exception):
    """Raised when query parsing fails"""
    pass


class QueryParser:
    """Parse and validate Pharmit query JSON format"""
    
    def __init__(self, query_path: Path):
        if not isinstance(query_path, Path):
            query_path = Path(query_path)
        
        if not query_path.exists():
            raise FileNotFoundError(f"Query file not found: {query_path}")
        
        self.query_path = query_path
        self.points = []
        self.num_points = 0
        
    def parse(self) -> bool:
        """Parse query JSON file and extract pharmacophore points"""
        try:
            with open(self.query_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            
            if 'points' not in data:
                raise QueryParseError("Query JSON must contain 'points' array")
            
            raw_points = data['points']
            
            for i, point in enumerate(raw_points):
                # Skip disabled points
                if not point.get('enabled', True):
                    logger.debug(f"Skipping disabled point {i}")
                    continue
                
                # Validate required fields
                if 'name' not in point:
                    raise QueryParseError(f"Point {i} missing 'name' field")
                
                if point['name'] not in PHARMACOPHORE_TYPES:
                    raise QueryParseError(f"Unknown pharmacophore type: {point['name']}")
                
                for coord in ['x', 'y', 'z']:
                    if coord not in point:
                        raise QueryParseError(f"Point {i} missing coordinate '{coord}'")
                
                # Create point entry
                parsed_point = {
                    'index': len(self.points),
                    'type': point['name'],
                    'type_id': PHARMACOPHORE_TYPES[point['name']],
                    'x': float(point['x']),
                    'y': float(point['y']),
                    'z': float(point['z']),
                    'coords': np.array([float(point['x']), float(point['y']), float(point['z'])]),
                    'radius': float(point.get('radius', 1.0)),
                    'requirement': point.get('requirement', 'Required')
                }
                
                # Check for H-bond vector constraint
                if point.get('vector_on', True) and 'svector' in point:
                    vec = point['svector']
                    vec_array = np.array([float(vec['x']), float(vec['y']), float(vec['z'])])
                    vec_norm = np.linalg.norm(vec_array)
                    if vec_norm > 0:
                        parsed_point['vector'] = vec_array / vec_norm
                
                self.points.append(parsed_point)
            
            self.num_points = len(self.points)
            
            if self.num_points < 3:
                raise QueryParseError(f"Need at least 3 points for triplet matching, got {self.num_points}")
            
            logger.info(f"Parsed {self.num_points} query pharmacophore points")
            return True
            
        except json.JSONDecodeError as e:
            raise QueryParseError(f"Invalid JSON: {e}")
        except (KeyError, ValueError, TypeError) as e:
            raise QueryParseError(f"Error parsing query: {e}")
    
    def get_points(self) -> List[Dict]:
        """Get parsed pharmacophore points"""
        return self.points


class TripletMatcher:
    """Find matching triplets in database using spatial indices"""
    
    def __init__(self, db_path: Path, distance_tolerance: float):
        if not isinstance(db_path, Path):
            db_path = Path(db_path)
        
        if not db_path.exists():
            raise FileNotFoundError(f"Database not found: {db_path}")
        
        self.db_path = db_path
        self.distance_tolerance = distance_tolerance
        self.conn = None
        
    def connect(self) -> None:
        """Connect to SQLite database"""
        try:
            self.conn = sqlite3.connect(str(self.db_path))
            self.conn.row_factory = sqlite3.Row
            
            # Verify schema
            cursor = self.conn.cursor()
            tables = cursor.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            ).fetchall()
            table_names = [t['name'] for t in tables]
            
            required_tables = ['molecules', 'features', 'triplets', 'pharma_types']
            missing = [t for t in required_tables if t not in table_names]
            if missing:
                raise DatabaseConnectionError(f"Missing required tables: {missing}")
            
            logger.info(f"Connected to database: {self.db_path}")
            
        except sqlite3.Error as e:
            raise DatabaseConnectionError(f"Database connection failed: {e}")
        
    def generate_query_triplets(self, points: List[Dict]) -> List[Dict]:
        """Generate all C(n,3) query triplets"""
        triplets = []
        
        for i, j, k in itertools.combinations(range(len(points)), 3):
            p1, p2, p3 = points[i], points[j], points[k]
            
            # Calculate pairwise distances
            d12 = np.linalg.norm(p1['coords'] - p2['coords'])
            d23 = np.linalg.norm(p2['coords'] - p3['coords'])
            d31 = np.linalg.norm(p3['coords'] - p1['coords'])
            
            # Canonical type ordering (sorted type IDs)
            types = sorted([
                (p1['type_id'], i),
                (p2['type_id'], j),
                (p3['type_id'], k)
            ])
            
            triplet = {
                'query_indices': [i, j, k],
                'type1_id': types[0][0],
                'type2_id': types[1][0],
                'type3_id': types[2][0],
                'd12': d12,
                'd23': d23,
                'd31': d31,
                'd12_min': d12 - self.distance_tolerance,
                'd12_max': d12 + self.distance_tolerance,
                'd23_min': d23 - self.distance_tolerance,
                'd23_max': d23 + self.distance_tolerance,
                'd31_min': d31 - self.distance_tolerance,
                'd31_max': d31 + self.distance_tolerance
            }
            
            triplets.append(triplet)
        
        logger.info(f"Generated {len(triplets)} query triplets")
        return triplets
    
    def rank_triplets(self, triplets: List[Dict]) -> List[Dict]:
        """Rank triplets by database frequency (rarest first for pruning)"""
        cursor = self.conn.cursor()
        
        for triplet in triplets:
            count = cursor.execute(
                "SELECT COUNT(*) as cnt FROM triplets WHERE type1_id=? AND type2_id=? AND type3_id=?",
                (triplet['type1_id'], triplet['type2_id'], triplet['type3_id'])
            ).fetchone()['cnt']
            
            triplet['db_count'] = count
            triplet['rarity'] = 1.0 / (count + 1.0)
        
        ranked = sorted(triplets, key=lambda x: x['rarity'], reverse=True)
        
        if ranked:
            logger.info(f"Rarest triplet has {ranked[0]['db_count']} database matches")
        
        return ranked
    
    def find_matches(self, query_triplets: List[Dict]) -> Dict[int, List[Dict]]:
        """Find database triplets matching query triplets within tolerance"""
        cursor = self.conn.cursor()
        matches_by_triplet = {}
        
        for idx, qt in enumerate(query_triplets):
            sql = """
            SELECT 
                t.triplet_id, t.mol_id,
                t.type1_id, t.type2_id, t.type3_id,
                t.dist12, t.dist23, t.dist31,
                t.point1_idx, t.point2_idx, t.point3_idx
            FROM triplets t
            WHERE 
                t.type1_id = ? AND
                t.type2_id = ? AND
                t.type3_id = ? AND
                t.dist12 >= ? AND t.dist12 <= ? AND
                t.dist23 >= ? AND t.dist23 <= ? AND
                t.dist31 >= ? AND t.dist31 <= ?
            """
            
            params = (
                qt['type1_id'], qt['type2_id'], qt['type3_id'],
                qt['d12_min'], qt['d12_max'],
                qt['d23_min'], qt['d23_max'],
                qt['d31_min'], qt['d31_max']
            )
            
            matches = [dict(row) for row in cursor.execute(sql, params)]
            matches_by_triplet[idx] = matches
            
            logger.debug(f"Query triplet {idx}: {len(matches)} database matches")
        
        total = sum(len(m) for m in matches_by_triplet.values())
        logger.info(f"Total triplet matches: {total}")
        
        return matches_by_triplet
    
    def group_by_molecule(self, matches_by_triplet: Dict[int, List[Dict]]) -> Dict[int, Dict]:
        """Group triplet matches by molecule ID"""
        mol_data = defaultdict(lambda: defaultdict(list))
        
        for triplet_idx, matches in matches_by_triplet.items():
            for match in matches:
                mol_id = match['mol_id']
                mol_data[mol_id][triplet_idx].append(match)
        
        logger.info(f"Matches span {len(mol_data)} unique molecules")
        return dict(mol_data)
    
    def close(self) -> None:
        """Close database connection"""
        if self.conn:
            self.conn.close()


class CorrespondenceFinder:
    """Recursive backtracking to find consistent pharmacophore correspondences"""
    
    def __init__(self, query_points: List[Dict], query_triplets: List[Dict], 
                 mol_matches: Dict, db_conn, check_vectors: bool = True):
        self.query_points = query_points
        self.query_triplets = query_triplets
        self.mol_matches = mol_matches
        self.db_conn = db_conn
        self.check_vectors = check_vectors
        self.results = []
        
    def search(self, max_results: int = 100) -> List[Dict]:
        """Find all consistent correspondences"""
        for mol_id, triplet_matches in self.mol_matches.items():
            # Start backtracking from first (rarest) triplet
            if 0 in triplet_matches:
                for db_triplet in triplet_matches[0]:
                    self._backtrack(
                        mol_id=mol_id,
                        triplet_idx=0,
                        correspondence={},
                        matched_db_points=set(),
                        triplet_matches=triplet_matches,
                        db_triplet=db_triplet
                    )
            
            if len(self.results) >= max_results:
                logger.info(f"Reached max results ({max_results}), stopping search")
                break
        
        logger.info(f"Found {len(self.results)} valid correspondences")
        return self.results
    
    def _backtrack(self, mol_id, triplet_idx, correspondence, 
                   matched_db_points, triplet_matches, db_triplet):
        """Recursive backtracking"""
        
        # Check if this triplet is consistent
        query_triplet = self.query_triplets[triplet_idx]
        query_indices = query_triplet['query_indices']
        db_indices = [db_triplet['point1_idx'], db_triplet['point2_idx'], db_triplet['point3_idx']]
        
        # Verify consistency
        new_mapping = {}
        for q_idx, db_idx in zip(query_indices, db_indices):
            if q_idx in correspondence:
                if correspondence[q_idx] != db_idx:
                    return  # Inconsistent
            else:
                if db_idx in matched_db_points:
                    return  # DB point already used
                new_mapping[q_idx] = db_idx
        
        # Check vectors if enabled
        if self.check_vectors and new_mapping:
            if not self._check_vectors(mol_id, db_triplet, query_indices, new_mapping):
                return
        
        # Update correspondence
        new_correspondence = correspondence.copy()
        new_correspondence.update(new_mapping)
        
        new_matched = matched_db_points.copy()
        new_matched.update(new_mapping.values())
        
        # Base case: all triplets matched
        if triplet_idx == len(self.query_triplets) - 1:
            if len(new_correspondence) == len(self.query_points):
                self.results.append({
                    'mol_id': mol_id,
                    'correspondence': new_correspondence,
                    'num_matched': len(new_correspondence)
                })
            return
        
        # Recursive case: try next triplet
        next_idx = triplet_idx + 1
        if next_idx in triplet_matches:
            for next_db_triplet in triplet_matches[next_idx]:
                self._backtrack(
                    mol_id, next_idx, new_correspondence,
                    new_matched, triplet_matches, next_db_triplet
                )
    
    def _check_vectors(self, mol_id, db_triplet, query_indices, new_mapping):
        """Check vector constraints for H-bonds"""
        cursor = self.db_conn.cursor()
        
        for q_idx in new_mapping:
            query_point = self.query_points[q_idx]
            
            # Skip if no vector constraint
            if 'vector' not in query_point:
                continue
            
            # Get DB feature vector
            db_point_idx = new_mapping[q_idx]
            
            # Query from features table
            row = cursor.execute("""
                SELECT vector_x, vector_y, vector_z
                FROM features
                WHERE mol_id = ? AND point_idx = ?
            """, (mol_id, db_point_idx)).fetchone()
            
            if not row:
                continue
            
            if row['vector_x'] is None:
                continue  # No vector stored
            
            db_vector = np.array([row['vector_x'], row['vector_y'], row['vector_z']])
            vec_norm = np.linalg.norm(db_vector)
            if vec_norm > 0:
                db_vector = db_vector / vec_norm
            else:
                continue
            
            # Calculate angle
            cos_angle = np.dot(query_point['vector'], db_vector)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            angle_deg = np.degrees(np.arccos(cos_angle))
            
            # Reject if angle > 60 degrees
            if angle_deg > 60:
                return False
        
        return True


class RMSDCalculator:
    """Calculate RMSD using Kabsch algorithm"""
    
    @staticmethod
    def kabsch_rmsd(query_coords: np.ndarray, db_coords: np.ndarray, 
                    weights: Optional[np.ndarray] = None) -> Tuple[float, np.ndarray, np.ndarray]:
        """
        Calculate weighted RMSD with optimal superposition using Kabsch algorithm
        
        Returns: (rmsd, rotation_matrix, translation_vector)
        """
        if weights is None:
            weights = np.ones(len(query_coords))
        
        # Normalize weights
        weights = weights / np.sum(weights)
        
        # Calculate weighted centroids
        query_centroid = np.average(query_coords, axis=0, weights=weights)
        db_centroid = np.average(db_coords, axis=0, weights=weights)
        
        # Center coordinates
        query_centered = query_coords - query_centroid
        db_centered = db_coords - db_centroid
        
        # Weighted covariance matrix
        H = np.dot((query_centered * weights[:, np.newaxis]).T, db_centered)
        
        # SVD
        U, S, Vt = np.linalg.svd(H)
        rotation = np.dot(Vt.T, U.T)
        
        # Ensure proper rotation (det = +1)
        if np.linalg.det(rotation) < 0:
            Vt[-1, :] *= -1
            rotation = np.dot(Vt.T, U.T)
        
        # Apply rotation
        db_rotated = np.dot(db_centered, rotation)
        
        # Calculate weighted RMSD
        diff = query_centered - db_rotated
        weighted_sq_dist = np.sum(diff**2, axis=1) * weights
        rmsd = np.sqrt(np.sum(weighted_sq_dist))
        
        translation = query_centroid - np.dot(db_centroid, rotation)
        
        return rmsd, rotation, translation
    
    @staticmethod
    def calculate_correspondence_rmsd(correspondence: Dict, query_points: List[Dict], 
                                     mol_id: int, db_conn) -> Optional[Dict]:
        """Calculate RMSD for a correspondence"""
        cursor = db_conn.cursor()
        
        # Extract coordinates
        query_coords = []
        db_coords = []
        weights = []
        
        for q_idx in sorted(correspondence.keys()):
            query_point = query_points[q_idx]
            db_point_idx = correspondence[q_idx]
            
            # Query coordinates
            query_coords.append(query_point['coords'])
            
            # Weight from radius (1/radius^2)
            weights.append(1.0 / (query_point['radius'] ** 2))
            
            # DB coordinates from features table
            row = cursor.execute("""
                SELECT x, y, z
                FROM features
                WHERE mol_id = ? AND point_idx = ?
            """, (mol_id, db_point_idx)).fetchone()
            
            if not row:
                return None
            
            db_coords.append(np.array([row['x'], row['y'], row['z']]))
        
        query_coords = np.array(query_coords)
        db_coords = np.array(db_coords)
        weights = np.array(weights)
        
        # Kabsch RMSD
        rmsd, rotation, translation = RMSDCalculator.kabsch_rmsd(
            query_coords, db_coords, weights
        )
        
        return {
            'rmsd': float(rmsd),
            'rotation': rotation.tolist(),
            'translation': translation.tolist(),
            'num_atoms': len(query_coords)
        }


class PharmacophoreSearch:
    """Main search engine orchestrator"""
    
    def __init__(self, query_path: Path, db_path: Path, 
                 distance_tolerance: float = 1.0,
                 rmsd_threshold: float = 2.0,
                 max_results: int = 20,
                 check_vectors: bool = True):
        self.query_path = query_path
        self.db_path = db_path
        self.distance_tolerance = distance_tolerance
        self.rmsd_threshold = rmsd_threshold
        self.max_results = max_results
        self.check_vectors = check_vectors
        self.results = []
        
    def search(self) -> List[Dict]:
        """Execute complete search pipeline"""
        logger.info("="*70)
        logger.info("PHARMACOPHORE DATABASE SEARCH")
        logger.info("="*70)
        logger.info(f"Query: {self.query_path}")
        logger.info(f"Database: {self.db_path}")
        logger.info(f"Distance tolerance: {self.distance_tolerance} Å")
        logger.info(f"RMSD threshold: {self.rmsd_threshold} Å")
        logger.info(f"Max results: {self.max_results}")
        logger.info("="*70)
        
        # Stage 1: Parse query
        logger.info("\nStage 1: Parsing query...")
        parser = QueryParser(self.query_path)
        if not parser.parse():
            logger.error("Failed to parse query")
            return []
        
        query_points = parser.get_points()
        
        # Stage 2: Triplet matching
        logger.info("\nStage 2: Triplet matching...")
        matcher = TripletMatcher(self.db_path, self.distance_tolerance)
        matcher.connect()
        
        query_triplets = matcher.generate_query_triplets(query_points)
        ranked_triplets = matcher.rank_triplets(query_triplets)
        matches_by_triplet = matcher.find_matches(ranked_triplets)
        mol_matches = matcher.group_by_molecule(matches_by_triplet)
        
        # Stage 3: Recursive backtracking
        logger.info("\nStage 3: Recursive backtracking...")
        finder = CorrespondenceFinder(
            query_points, ranked_triplets, mol_matches,
            matcher.conn, self.check_vectors
        )
        correspondences = finder.search(self.max_results * 10)  # Get more candidates
        
        # Stage 4: RMSD calculation and ranking
        logger.info("\nStage 4: RMSD calculation and ranking...")
        scored_results = []
        
        for corr in correspondences:
            rmsd_result = RMSDCalculator.calculate_correspondence_rmsd(
                corr['correspondence'], query_points,
                corr['mol_id'], matcher.conn
            )
            
            if rmsd_result and rmsd_result['rmsd'] <= self.rmsd_threshold:
                # Get molecule info
                mol_info = matcher.conn.execute("""
                    SELECT pdb_id, ligand_name, chain, num_features
                    FROM molecules WHERE mol_id = ?
                """, (corr['mol_id'],)).fetchone()
                
                result = {
                    'pdb_id': mol_info['pdb_id'],
                    'ligand_name': mol_info['ligand_name'],
                    'chain': mol_info['chain'] or '',
                    'rmsd': rmsd_result['rmsd'],
                    'num_matched_features': corr['num_matched'],
                    'num_query_features': len(query_points),
                    'num_db_features': mol_info['num_features'],
                    'match_percentage': 100.0 * corr['num_matched'] / len(query_points)
                }
                
                scored_results.append(result)
        
        # Sort by RMSD (best first) and limit
        scored_results.sort(key=lambda x: x['rmsd'])
        self.results = scored_results[:self.max_results]
        
        matcher.close()
        
        logger.info(f"\nFound {len(self.results)} matches (RMSD ≤ {self.rmsd_threshold} Å)")
        logger.info("="*70)
        
        return self.results
    
    def save_results(self, output_csv: Path, log_file: Path):
        """Save results to CSV and log file"""
        if not self.results:
            logger.warning("No results to save")
            return
        
        # Create DataFrame
        df = pd.DataFrame(self.results)
        df.insert(0, 'rank', range(1, len(df) + 1))
        
        # Save CSV
        df.to_csv(output_csv, index=False, float_format='%.4f')
        logger.info(f"Results saved to: {output_csv}")
        
        # Save detailed log
        with open(log_file, 'w', encoding='utf-8') as f:
            f.write("PHARMACOPHORE DATABASE SEARCH RESULTS\n")
            f.write("="*70 + "\n\n")
            f.write(f"Query: {self.query_path}\n")
            f.write(f"Database: {self.db_path}\n")
            f.write(f"Search date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            f.write(f"Parameters:\n")
            f.write(f"  Distance tolerance: {self.distance_tolerance} Å\n")
            f.write(f"  RMSD threshold: {self.rmsd_threshold} Å\n")
            f.write(f"  Max results: {self.max_results}\n")
            f.write(f"  Vector constraints: {self.check_vectors}\n\n")
            f.write("="*70 + "\n\n")
            f.write(f"Total matches found: {len(self.results)}\n\n")
            f.write("Top matches:\n")
            f.write("-"*70 + "\n")
            
            for i, result in enumerate(self.results, 1):
                f.write(f"\nRank {i}:\n")
                f.write(f"  PDB ID: {result['pdb_id']}\n")
                f.write(f"  Ligand: {result['ligand_name']}\n")
                f.write(f"  Chain: {result['chain']}\n")
                f.write(f"  RMSD: {result['rmsd']:.4f} Å\n")
                f.write(f"  Matched features: {result['num_matched_features']}/{result['num_query_features']}\n")
                f.write(f"  Match percentage: {result['match_percentage']:.1f}%\n")
        
        logger.info(f"Detailed log saved to: {log_file}")


def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(
        description='Search pharmacophore database for matching molecules',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  python search_pharmacophore_database.py --query myquery.json --database Database/pharmacophores.db
  
  python search_pharmacophore_database.py \\
    --query query.json \\
    --database pharmacophores.db \\
    --output results.csv \\
    --log search.log \\
    --distance-tolerance 1.0 \\
    --rmsd-threshold 2.0 \\
    --max-results 20
        '''
    )
    
    parser.add_argument('--query', type=str, required=True,
                       help='Query pharmacophore JSON file')
    parser.add_argument('--database', type=str, required=True,
                       help='Pharmacophore database file (.db)')
    parser.add_argument('--output', type=str, default='search_results.csv',
                       help='Output CSV file (default: search_results.csv)')
    parser.add_argument('--log', type=str, default='search.log',
                       help='Log file (default: search.log)')
    parser.add_argument('--distance-tolerance', type=float, default=1.0,
                       help='Distance tolerance in Angstroms (default: 1.0)')
    parser.add_argument('--rmsd-threshold', type=float, default=2.0,
                       help='RMSD threshold in Angstroms (default: 2.0)')
    parser.add_argument('--max-results', type=int, default=20,
                       help='Maximum number of results (default: 20)')
    parser.add_argument('--ignore-vectors', action='store_true',
                       help='Ignore H-bond vector constraints')
    parser.add_argument('--verbose', action='store_true',
                       help='Verbose logging')
    
    args = parser.parse_args()
    
    # Setup logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler('search_debug.log'),
            logging.StreamHandler(sys.stdout)
        ]
    )
    
    # Validate inputs
    query_path = Path(args.query)
    db_path = Path(args.database)
    
    if not query_path.exists():
        logger.error(f"Query file not found: {query_path}")
        return 1
    
    if not db_path.exists():
        logger.error(f"Database file not found: {db_path}")
        return 1
    
    try:
        # Execute search
        search = PharmacophoreSearch(
            query_path=query_path,
            db_path=db_path,
            distance_tolerance=args.distance_tolerance,
            rmsd_threshold=args.rmsd_threshold,
            max_results=args.max_results,
            check_vectors=not args.ignore_vectors
        )
        
        results = search.search()
        
        if results:
            output_csv = Path(args.output)
            log_file = Path(args.log)
            search.save_results(output_csv, log_file)
            logger.info(f"\n✓ Search complete! Found {len(results)} matches")
            return 0
        else:
            logger.warning("\n✗ No matches found")
            return 0  # Not an error, just no results
            
    except (QueryParseError, DatabaseConnectionError) as e:
        logger.error(f"\n✗ Error: {e}")
        return 1
    except Exception as e:
        logger.error(f"\n✗ Unexpected error: {e}", exc_info=True)
        return 1


if __name__ == '__main__':
    sys.exit(main())
