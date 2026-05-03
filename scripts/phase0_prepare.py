import os
import urllib.request
from Bio.PDB import PDBParser, PDBIO, Structure, Model, Chain
import warnings
from Bio import BiopythonWarning
from utils.logger import console

warnings.simplefilter('ignore', BiopythonWarning)

def run_phase0(input_dir: str):
    console.rule("[bold cyan]Phase 0: Target Data Preparation[/bold cyan]")
    output_pdb = os.path.join(input_dir, "rsvf_3_motifs.pdb")
    
    if os.path.exists(output_pdb):
        console.print(f"[success]✔[/success] Target motifs already exist at {output_pdb}. Skipping download.")
        return output_pdb
        
    console.print("[info]Downloading PDB 4JHW (RSVF Prefusion)...[/info]")
    url = "https://files.rcsb.org/download/4JHW.pdb"
    temp_pdb = os.path.join(input_dir, "4JHW.pdb")
    if not os.path.exists(temp_pdb):
        urllib.request.urlretrieve(url, temp_pdb)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("4JHW", temp_pdb)
    chain_f = structure[0]['F']

    motifs = {
        'A': {'name': 'Site_II', 'start': 254, 'end': 277}, 
        'B': {'name': 'Site_IV', 'start': 422, 'end': 438}, 
        'C': {'name': 'Site_V',  'start': 163, 'end': 181}, 
    }

    new_struct = Structure.Structure("RSVF_Motifs")
    new_model = Model.Model(0)
    new_struct.add(new_model)

    for new_chain_id, info in motifs.items():
        new_chain = Chain.Chain(new_chain_id)
        console.print(f"  [info]► Extracting {info['name']} to Chain {new_chain_id}...[/info]")
        for res_id in range(info['start'], info['end'] + 1):
            if (' ', res_id, ' ') in chain_f:
                res = chain_f[(' ', res_id, ' ')].copy()
                res.id = (' ', res_id - info['start'] + 1, ' ') 
                new_chain.add(res)
        new_model.add(new_chain)

    io = PDBIO()
    io.set_structure(new_struct)
    io.save(output_pdb)
    console.print(f"[success]✔ Motifs successfully packaged to {output_pdb}[/success]")
    return output_pdb
