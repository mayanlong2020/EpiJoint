import os
from Bio.PDB import PDBParser
import numpy as np

parser = PDBParser(QUIET=True)
native_pdb = os.path.expanduser("~/Development/Codebase/BakerLab/inputs/rsvf_3_motifs.pdb")
struct = parser.get_structure("native", native_pdb)
atoms = [atom for atom in struct.get_atoms() if atom.get_name() == 'CA']

# Group atoms by chain
chains = {}
for atom in atoms:
    chain_id = atom.get_parent().get_parent().id
    if chain_id not in chains:
        chains[chain_id] = []
    chains[chain_id].append(atom.get_coord())

all_coords = np.array([a.get_coord() for a in atoms])
center_all = np.mean(all_coords, axis=0)

print(f"Center of all motifs: {center_all}")
for chain_id, coords in chains.items():
    coords = np.array(coords)
    center_motif = np.mean(coords, axis=0)
    vec = center_motif - center_all
    print(f"Motif {chain_id}: Center {center_motif}, Vector from center {vec}, Norm {np.linalg.norm(vec):.2f}")
