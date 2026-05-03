import pickle
import os

run_dir = os.path.expanduser("~/Development/Codebase/BakerLab/runs/run_20260503_012538")
backbone_dir = os.path.join(run_dir, "01_backbones")
trb_files = [f for f in os.listdir(backbone_dir) if f.endswith(".trb")]

if trb_files:
    trb_path = os.path.join(backbone_dir, trb_files[0])
    with open(trb_path, 'rb') as f:
        trb = pickle.load(f)
    print("Keys in trb:", trb.keys())
    print("con_ref_pdb_idx:", trb.get('con_ref_pdb_idx'))
    print("con_hal_pdb_idx:", trb.get('con_hal_pdb_idx'))
    # check for motif ranges if available
    if 'con_ref_idx0' in trb:
        print("con_ref_idx0:", trb['con_ref_idx0'])
    if 'con_hal_idx0' in trb:
        print("con_hal_idx0:", trb['con_hal_idx0'])
else:
    print("No trb files found.")
