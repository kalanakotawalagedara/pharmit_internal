#!/usr/bin/env python3
"""
Interaction Pharmacophore Generator
Generates protein-ligand interaction pharmacophore models from PDB structures
Based on Pharmit's interaction detection algorithms
"""

import os
import sys
import json
import logging
import traceback
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import pandas as pd
import requests
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
import numpy as np


# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('pharmacophore_generation.log'),
        logging.StreamHandler(sys.stdout)
    ]
)
logger = logging.getLogger(__name__)


# Pharmit SMARTS patterns from pharmarec.cpp (exact C++ definitions)
PHARMACOPHORE_DEFINITIONS = {
    'Aromatic': {
        'smarts': ['a1aaaaa1', 'a1aaaa1'],
        'radius': 1.0,
        'has_vector': False
    },
    'HydrogenDonor': {
        'smarts': [
            '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]',
            '[#8!H0&!$([OH][C,S,P]=O)]',
            '[#16!H0]'
        ],
        'radius': 1.0,
        'has_vector': True
    },
    'HydrogenAcceptor': {
        'smarts': [
            '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]',
            '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]'
        ],
        'radius': 1.0,
        'has_vector': True
    },
    'PositiveIon': {
        'smarts': [
            '[+,+2,+3,+4]',
            '[$(CC)](=N)N',  # amidine
            '[$(C(N)(N)=N)]',  # guanidine
            '[$(n1cc[nH]c1)]'
        ],
        'radius': 1.0,
        'has_vector': False
    },
    'NegativeIon': {
        'smarts': [
            '[-,-2,-3,-4]',
            'C(=O)[O-,OH,OX1]',
            '[$([S,P](=O)[O-,OH,OX1])]',
            'c1[nH1]nnn1',
            'c1nn[nH1]n1',
            'C(=O)N[OH1,O-,OX1]',
            'C(=O)N[OH1,O-]',
            'CO(=N[OH1,O-])',
            '[$(N-[SX4](=O)(=O)[CX4](F)(F)F)]'  # trifluoromethyl sulfonamide
        ],
        'radius': 1.0,
        'has_vector': False
    },
    'Hydrophobic': {
        'smarts': [
            'a1aaaaa1',
            'a1aaaa1',
            '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]',
            '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[C&r3]1~[C&r3]~[C&r3]1',
            '[C&r4]1~[C&r4]~[C&r4]~[C&r4]1',
            '[C&r5]1~[C&r5]~[C&r5]~[C&r5]~[C&r5]1',
            '[C&r6]1~[C&r6]~[C&r6]~[C&r6]~[C&r6]~[C&r6]1',
            '[C&r7]1~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]~[C&r7]1',
            '[C&r8]1~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]~[C&r8]1',
            '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]',
            '[$([CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]',
            '[$([S]~[#6])&!$(S~[!#6])]'
        ],
        'radius': 1.0,
        'has_vector': False
    }
}

# Protein SMARTS patterns (optimized for receptor feature detection)
PROTEIN_PHARMACOPHORE_DEFINITIONS = {
    'Aromatic': {
        'smarts': ['a1aaaaa1', 'a1aaaa1', '[+,+2,+3,+4]', '[$(C(N)(N)=N)]', '[$(n1cc[nH]c1)]'],
        'radius': 1.0
    },
    'HydrogenDonor': {
        'smarts': [
            '[#7!H0&!$(N-[SX4](=O)(=O)[CX4](F)(F)F)]',
            '[#8!H0&!$([OH][C,S,P]=O)]',
            '[#16!H0]'
        ],
        'radius': 1.0
    },
    'HydrogenAcceptor': {
        'smarts': [
            '[#7&!$([nX3])&!$([NX3]-*=[!#6])&!$([NX3]-[a])&!$([NX4])&!$(N=C([C,N])N)]',
            '[$([O])&!$([OX2](C)C=O)&!$(*(~a)~a)]'
        ],
        'radius': 1.0
    },
    'PositiveIon': {
        'smarts': ['[+,+2,+3,+4]', '[$(C(N)(N)=N)]', '[$(n1cc[nH]c1)]'],
        'radius': 1.0
    },
    'NegativeIon': {
        'smarts': ['[-,-2,-3,-4]', 'C(=O)[O-,OH,OX1]'],
        'radius': 1.0
    },
    'Hydrophobic': {
        'smarts': [
            'a1aaaaa1',
            'a1aaaa1',
            '[$([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(**[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]',
            '[$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])&!$(*([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I])]([CH3X4,CH2X3,CH1X2,F,Cl,Br,I])[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[CH2X4,CH1X3,CH0X2]~[CH3X4,CH2X3,CH1X2,F,Cl,Br,I]',
            '[$([CH2X4,CH1X3,CH0X2]~[$([!#1]);!$([CH2X4,CH1X3,CH0X2])])]~[CH2X4,CH1X3,CH0X2]~[CH2X4,CH1X3,CH0X2]',
            '[$([S]~[#6])&!$(S~[!#6])]'
        ],
        'radius': 1.0
    }
}

# Interaction rules (from pharmaInteractions map in C++)
# Format: 'LigandType': (complementReceptorType, maxDistance, minMatches)
INTERACTION_RULES = {
    'Aromatic': ('Aromatic', 5.0, 1),
    'HydrogenDonor': ('HydrogenAcceptor', 4.0, 1),
    'HydrogenAcceptor': ('HydrogenDonor', 4.0, 1),
    'PositiveIon': ('NegativeIon', 5.0, 1),
    'NegativeIon': ('PositiveIon', 5.0, 1),
    'Hydrophobic': ('Hydrophobic', 6.0, 3)
}


class PDBParser:
    """Parse PDB files and extract ligand/receptor structures"""
    
    @staticmethod
    def download_pdb(pdb_id: str) -> Optional[str]:
        """Download PDB file from RCSB"""
        url = f"https://files.rcsb.org/view/{pdb_id.upper()}.pdb"
        try:
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            return response.text
        except Exception as e:
            logger.error(f"Failed to download {pdb_id}: {e}")
            return None
    
    @staticmethod
    def parse_pdb_content(pdb_content: str, ligand_name: str) -> Tuple[Optional[str], Optional[str], Optional[str], Optional[str]]:
        """
        Parse PDB content and extract ligand and receptor
        Returns: (ligand_pdb, receptor_pdb, ligand_chain, water_pdb)
        """
        lines = pdb_content.split('\n')
        ligand_lines = []
        receptor_lines = []
        water_lines = []
        ligand_chain = None
        ligand_coords = []
        
        # First pass: extract ligand
        for line in lines:
            if not (line.startswith('HETATM') or line.startswith('ATOM')):
                continue
                
            resname = line[17:20].strip()
            chain = line[21:22]
            
            if resname == ligand_name:
                if ligand_chain is None:
                    ligand_chain = chain
                if chain == ligand_chain:
                    ligand_lines.append(line)
                    # Extract coordinates
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                        ligand_coords.append([x, y, z])
                    except ValueError:
                        continue
            elif line.startswith('ATOM'):
                receptor_lines.append(line)
            elif resname in ['HOH', 'WAT']:
                water_lines.append(line)
        
        if not ligand_lines:
            return None, None, None, None
        
        # Calculate ligand bounding box (with 4Å extension for water detection)
        if ligand_coords:
            coords_array = np.array(ligand_coords)
            min_coords = coords_array.min(axis=0) - 4.0
            max_coords = coords_array.max(axis=0) + 4.0
            
            # Filter binding site waters
            binding_waters = []
            for line in water_lines:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    
                    if (min_coords[0] <= x <= max_coords[0] and
                        min_coords[1] <= y <= max_coords[1] and
                        min_coords[2] <= z <= max_coords[2]):
                        # Check if close to any ligand atom
                        for lc in ligand_coords:
                            dist_sq = sum((lc[i] - [x,y,z][i])**2 for i in range(3))
                            if dist_sq <= 16.0:  # 4Å squared
                                binding_waters.append(line)
                                break
                except ValueError:
                    continue
        else:
            binding_waters = []
        
        ligand_pdb = '\n'.join(ligand_lines)
        receptor_pdb = '\n'.join(receptor_lines)
        water_pdb = '\n'.join(binding_waters)
        
        return ligand_pdb, receptor_pdb, ligand_chain, water_pdb


class PharmacophoreDetector:
    """Detect pharmacophore features from molecules"""
    
    def __init__(self):
        # Compile SMARTS patterns
        self.ligand_patterns = {}
        self.protein_patterns = {}
        
        for feature_type, definition in PHARMACOPHORE_DEFINITIONS.items():
            self.ligand_patterns[feature_type] = [
                Chem.MolFromSmarts(smarts) for smarts in definition['smarts']
            ]
        
        for feature_type, definition in PROTEIN_PHARMACOPHORE_DEFINITIONS.items():
            self.protein_patterns[feature_type] = [
                Chem.MolFromSmarts(smarts) for smarts in definition['smarts']
            ]
    
    def detect_features(self, mol: Chem.Mol, is_protein: bool = False) -> List[Dict]:
        """Detect pharmacophore features in molecule"""
        if mol is None:
            return []
        
        features = []
        patterns = self.protein_patterns if is_protein else self.ligand_patterns
        conf = mol.GetConformer()
        
        for feature_type, smarts_list in patterns.items():
            matched_atoms = set()
            
            for smarts_mol in smarts_list:
                if smarts_mol is None:
                    continue
                    
                matches = mol.GetSubstructMatches(smarts_mol)
                for match in matches:
                    # Calculate centroid of matched atoms
                    if not match:
                        continue
                    
                    # Avoid duplicate features at same location
                    atom_tuple = tuple(sorted(match))
                    if atom_tuple in matched_atoms:
                        continue
                    matched_atoms.add(atom_tuple)
                    
                    coords = [conf.GetAtomPosition(idx) for idx in match]
                    centroid = np.mean([[c.x, c.y, c.z] for c in coords], axis=0)
                    
                    feature = {
                        'type': feature_type,
                        'x': float(centroid[0]),
                        'y': float(centroid[1]),
                        'z': float(centroid[2]),
                        'radius': PHARMACOPHORE_DEFINITIONS[feature_type]['radius'],
                        'atom_indices': list(match)
                    }
                    features.append(feature)
        
        return features


class InteractionPharmacophoreGenerator:
    """Generate interaction pharmacophore models"""
    
    def __init__(self):
        self.detector = PharmacophoreDetector()
    
    def calculate_distance(self, point1: Dict, point2: Dict) -> float:
        """Calculate Euclidean distance between two points"""
        dx = point1['x'] - point2['x']
        dy = point1['y'] - point2['y']
        dz = point1['z'] - point2['z']
        return np.sqrt(dx*dx + dy*dy + dz*dz)
    
    def calculate_unit_vector(self, from_point: Dict, to_point: Dict) -> Dict:
        """Calculate unit vector from one point to another"""
        dx = to_point['x'] - from_point['x']
        dy = to_point['y'] - from_point['y']
        dz = to_point['z'] - from_point['z']
        length = np.sqrt(dx*dx + dy*dy + dz*dz)
        
        if length < 0.0001:
            return {'x': 0.0, 'y': 0.0, 'z': 0.0}
        
        return {
            'x': float(dx / length),
            'y': float(dy / length),
            'z': float(dz / length)
        }
    
    def filter_interaction_features(
        self,
        ligand_features: List[Dict],
        receptor_features: List[Dict]
    ) -> Tuple[List[Dict], List[Dict]]:
        """
        Filter ligand features based on receptor interactions
        Returns: (validated_features, screened_out_features)
        """
        # Group receptor features by type
        receptor_by_type = {}
        for feature in receptor_features:
            ftype = feature['type']
            if ftype not in receptor_by_type:
                receptor_by_type[ftype] = []
            receptor_by_type[ftype].append(feature)
        
        validated_features = []
        screened_out = []
        
        for lig_feature in ligand_features:
            lig_type = lig_feature['type']
            
            # Check if this type has interaction rules
            if lig_type not in INTERACTION_RULES:
                screened_out.append(lig_feature)
                continue
            
            complement_type, max_dist, min_matches = INTERACTION_RULES[lig_type]
            
            # Find complementary receptor features
            if complement_type not in receptor_by_type:
                screened_out.append(lig_feature)
                continue
            
            complement_features = receptor_by_type[complement_type]
            match_count = 0
            closest_receptor = None
            min_distance = float('inf')
            
            for rec_feature in complement_features:
                dist = self.calculate_distance(lig_feature, rec_feature)
                if dist <= max_dist:
                    match_count += 1
                    if dist < min_distance:
                        min_distance = dist
                        closest_receptor = rec_feature
            
            # Validate based on minimum matches
            if match_count >= min_matches:
                validated_feature = lig_feature.copy()
                
                # Add directional vector for H-bond donors and acceptors
                if PHARMACOPHORE_DEFINITIONS[lig_type]['has_vector'] and closest_receptor:
                    vector = self.calculate_unit_vector(lig_feature, closest_receptor)
                    validated_feature['vector'] = vector
                    # Store interaction distance for deduplication
                    validated_feature['interaction_distance'] = min_distance
                
                validated_features.append(validated_feature)
            else:
                screened_out.append(lig_feature)
        
        return validated_features, screened_out
    
    def deduplicate_overlapping_features(self, features: List[Dict]) -> List[Dict]:
        """
        Remove duplicate features at the same position.
        For H-bond features (HydrogenDonor/HydrogenAcceptor) at the same coordinates,
        keep only the feature with the shortest interaction distance.
        """
        if not features:
            return features
        
        # Group features by rounded coordinates (0.01Å precision)
        position_groups = {}
        for feature in features:
            # Round to 2 decimal places
            pos_key = (
                round(feature['x'], 2),
                round(feature['y'], 2),
                round(feature['z'], 2)
            )
            if pos_key not in position_groups:
                position_groups[pos_key] = []
            position_groups[pos_key].append(feature)
        
        # Process each position group
        deduplicated = []
        for pos_key, group in position_groups.items():
            if len(group) == 1:
                # No duplicates at this position
                deduplicated.append(group[0])
            else:
                # Check if we have H-bond duplicates
                hbond_types = {'HydrogenDonor', 'HydrogenAcceptor'}
                hbond_features = [f for f in group if f['type'] in hbond_types]
                other_features = [f for f in group if f['type'] not in hbond_types]
                
                if len(hbond_features) > 1:
                    # Multiple H-bond features at same position - keep the one with shortest distance
                    best_hbond = min(hbond_features, key=lambda f: f.get('interaction_distance', float('inf')))
                    deduplicated.append(best_hbond)
                    deduplicated.extend(other_features)
                else:
                    # No H-bond duplicates, keep all features
                    deduplicated.extend(group)
        
        return deduplicated
    
    def generate_pharmacophore(
        self,
        ligand_pdb: str,
        receptor_pdb: str,
        water_pdb: Optional[str] = None
    ) -> Optional[Dict]:
        """Generate interaction pharmacophore from PDB strings"""
        try:
            # Parse ligand
            ligand_mol = Chem.MolFromPDBBlock(ligand_pdb, sanitize=True, removeHs=False)
            if ligand_mol is None:
                logger.error("Failed to parse ligand PDB")
                return None
            
            # Parse receptor (include waters as receptor atoms)
            receptor_content = receptor_pdb
            if water_pdb:
                receptor_content += '\n' + water_pdb
            
            receptor_mol = Chem.MolFromPDBBlock(receptor_content, sanitize=True, removeHs=False)
            if receptor_mol is None:
                logger.error("Failed to parse receptor PDB")
                return None
            
            # Detect features
            ligand_features = self.detector.detect_features(ligand_mol, is_protein=False)
            receptor_features = self.detector.detect_features(receptor_mol, is_protein=True)
            
            logger.info(f"Detected {len(ligand_features)} ligand features, {len(receptor_features)} receptor features")
            
            # Filter based on interactions
            validated_features, screened_out = self.filter_interaction_features(
                ligand_features, receptor_features
            )
            
            logger.info(f"Validated {len(validated_features)} interaction features (before deduplication)")
            
            # Remove duplicate H-bond features at same coordinates
            validated_features = self.deduplicate_overlapping_features(validated_features)
            
            logger.info(f"Final {len(validated_features)} unique interaction features")
            
            # Clean up temporary fields used for deduplication
            for feature in validated_features:
                feature.pop('interaction_distance', None)
            
            return {
                'validated_features': validated_features,
                'screened_features': screened_out,
                'ligand_atom_count': ligand_mol.GetNumAtoms(),
                'receptor_atom_count': receptor_mol.GetNumAtoms()
            }
            
        except Exception as e:
            logger.error(f"Error generating pharmacophore: {e}")
            logger.error(traceback.format_exc())
            return None


def process_entry(pdb_id: str, ligand_name: str, output_dir: Path) -> bool:
    """Process a single PDB entry"""
    try:
        logger.info(f"Processing {pdb_id}_{ligand_name}")
        
        # Download PDB
        pdb_content = PDBParser.download_pdb(pdb_id)
        if not pdb_content:
            logger.error(f"Failed to download {pdb_id}")
            return False
        
        # Parse PDB
        ligand_pdb, receptor_pdb, ligand_chain, water_pdb = PDBParser.parse_pdb_content(
            pdb_content, ligand_name
        )
        
        if not ligand_pdb or not receptor_pdb:
            logger.error(f"Failed to extract ligand {ligand_name} from {pdb_id}")
            return False
        
        logger.info(f"Extracted ligand from chain {ligand_chain}, {len(water_pdb.split(chr(10))) if water_pdb else 0} binding site waters")
        
        # Generate pharmacophore
        generator = InteractionPharmacophoreGenerator()
        result = generator.generate_pharmacophore(ligand_pdb, receptor_pdb, water_pdb)
        
        if not result:
            logger.error(f"Failed to generate pharmacophore for {pdb_id}_{ligand_name}")
            return False
        
        # Prepare output JSON
        output_data = {
            'pdb_id': pdb_id.upper(),
            'ligand_name': ligand_name,
            'ligand_chain': ligand_chain,
            'processing_date': datetime.utcnow().isoformat() + 'Z',
            'num_features': len(result['validated_features']),
            'features': result['validated_features'],
            'metadata': {
                'ligand_atoms': result['ligand_atom_count'],
                'receptor_atoms': result['receptor_atom_count'],
                'screened_out_features': len(result['screened_features'])
            }
        }
        
        # Save JSON
        output_file = output_dir / f"{pdb_id.upper()}_{ligand_name}.json"
        with open(output_file, 'w') as f:
            json.dump(output_data, f, indent=2)
        
        logger.info(f"Successfully generated {output_file}")
        return True
        
    except Exception as e:
        logger.error(f"Error processing {pdb_id}_{ligand_name}: {e}")
        logger.error(traceback.format_exc())
        return False


def main():
    """Main execution function"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate interaction pharmacophore models from PDB structures')
    parser.add_argument('--csv', type=str, required=True, help='Input CSV file with PDB_ID and Ligand_Name columns')
    parser.add_argument('--output', type=str, default='Output', help='Output directory for JSON files')
    parser.add_argument('--batch-start', type=int, default=0, help='Batch start index (for parallel processing)')
    parser.add_argument('--batch-size', type=int, default=None, help='Batch size (process subset of CSV)')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Read CSV
    try:
        df = pd.read_csv(args.csv)
        logger.info(f"Loaded {len(df)} entries from {args.csv}")
    except Exception as e:
        logger.error(f"Failed to read CSV: {e}")
        return 1
    
    # Validate columns
    if 'PDB_ID' not in df.columns or 'Ligand_Name' not in df.columns:
        logger.error("CSV must contain 'PDB_ID' and 'Ligand_Name' columns")
        return 1
    
    # Select batch if specified
    if args.batch_size:
        start_idx = args.batch_start
        end_idx = min(start_idx + args.batch_size, len(df))
        df = df.iloc[start_idx:end_idx]
        logger.info(f"Processing batch: rows {start_idx} to {end_idx-1}")
    
    # Process entries
    success_count = 0
    fail_count = 0
    error_log = []
    
    for idx, row in df.iterrows():
        pdb_id = str(row['PDB_ID']).strip()
        ligand_name = str(row['Ligand_Name']).strip()
        
        # Check if already processed (incremental mode)
        output_file = output_dir / f"{pdb_id.upper()}_{ligand_name}.json"
        if output_file.exists():
            logger.info(f"Skipping {pdb_id}_{ligand_name} (already processed)")
            success_count += 1
            continue
        
        success = process_entry(pdb_id, ligand_name, output_dir)
        
        if success:
            success_count += 1
        else:
            fail_count += 1
            error_log.append({
                'pdb_id': pdb_id,
                'ligand_name': ligand_name,
                'index': int(idx)
            })
    
    # Save error log
    if error_log:
        error_file = output_dir / 'error_log.json'
        with open(error_file, 'w') as f:
            json.dump(error_log, f, indent=2)
        logger.info(f"Error log saved to {error_file}")
    
    # Summary
    logger.info(f"\n{'='*50}")
    logger.info(f"Processing complete!")
    logger.info(f"Success: {success_count}")
    logger.info(f"Failed: {fail_count}")
    logger.info(f"Total: {success_count + fail_count}")
    logger.info(f"{'='*50}\n")
    
    return 0 if fail_count == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
