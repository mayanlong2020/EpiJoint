import os, random, itertools, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils.logger import console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn

os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"

RFJOINT_DIR = os.path.expanduser("~/Development/BioDesign/RFjoint2")

def run_phase1(run_dir: str, motif_pdb_path: str, num_designs: int = 20, max_workers: int = 14):
    console.rule("[bold cyan]Phase 1: Backbone Generation (RFjoint2)[/bold cyan]")
    
    out_dir = os.path.join(run_dir, "01_backbones")
    os.makedirs(out_dir, exist_ok=True)
    
    motifs = {'A': 24, 'B': 17, 'C': 19}
    permutations = list(itertools.permutations(['A', 'B', 'C']))
    
    console.print(f"[info]Engine:[ /info] 14-core Apple M4 Max | Concurrency: {max_workers}")
    console.print(f"[info]Target:[ /info] Generating {num_designs} Multi-motif Backbones\n")

    def worker(task_id):
        out_prefix = os.path.join(out_dir, f"design_{task_id}")
        # 💡 断点续传逻辑：如果 PDB 已存在，直接跳过
        if os.path.exists(f"{out_prefix}_0.pdb"):
            return True, task_id, "Skipped", ""

        order = random.choice(permutations)
        L1, L2, L3, L4 =[random.randint(10, 20) for _ in range(4)]
        contigs = f"{L1},{order[0]}1-{motifs[order[0]]},{L2},{order[1]}1-{motifs[order[1]]},{L3},{order[2]}1-{motifs[order[2]]},{L4}"
        abs_pdb_path = os.path.abspath(motif_pdb_path)
        
        cmd =[
            "conda", "run", "-n", "rfjoint2_m4",
            "python", "inpaint.py", 
            "--pdb", abs_pdb_path, 
            "--contigs", contigs, 
            "--out", out_prefix
        ]
        result = subprocess.run(cmd, cwd=RFJOINT_DIR, capture_output=True, text=True)
        
        if os.path.exists(f"{out_prefix}_0.pdb"):
            return True, task_id, contigs, ""
        else:
            return False, task_id, contigs, result.stderr

    success_count = 0
    with Progress(SpinnerColumn(), TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), TimeElapsedColumn(), console=console) as progress:
        task_pb = progress.add_task("[cyan]Forging backbones...", total=num_designs)
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(worker, i): i for i in range(num_designs)}
            for future in as_completed(futures):
                success, tid, contigs, err = future.result()
                if success: success_count += 1
                else: progress.console.print(f"[error]Task {tid} failed![/error] {err[-100:]}")
                progress.advance(task_pb)

    if success_count > 0:
        console.print(f"\n[success]✔ Phase 1 Completed! Forged {success_count}/{num_designs} backbones.[/success]")
    else:
        console.print(f"\n[error]✖ Phase 1 Failed.[/error]")
