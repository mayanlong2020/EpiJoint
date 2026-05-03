import os, time
from datetime import datetime
from utils.logger import console
from scripts.phase0_prepare import run_phase0
from scripts.phase1_generate import run_phase1

def main():
    console.print("\n[bold highlight]=================================================[/bold highlight]")
    console.print("[bold highlight]   🧬 BakerLab Pipeline: Multi-Epitope Forge   [/bold highlight]")
    console.print("[bold highlight]=================================================[/bold highlight]\n")
    
    base_dir = os.path.abspath(os.path.dirname(__file__))
    inputs_dir = os.path.join(base_dir, "inputs")
    runs_dir = os.path.join(base_dir, "runs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(runs_dir, exist_ok=True)
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    current_run_dir = os.path.join(runs_dir, f"run_{timestamp}")
    os.makedirs(current_run_dir, exist_ok=True)
    
    console.print(f"[info]Initialized Session:[/info] {current_run_dir}\n")
    start_time = time.time()
    
    motif_pdb = run_phase0(inputs_dir)
    run_phase1(current_run_dir, motif_pdb, num_designs=20, max_workers=14)
    
    total_time = time.time() - start_time
    console.rule(f"[bold highlight]Pipeline Paused. Time Elapsed: {total_time:.1f}s[/bold highlight]")
    console.print(f"\n[info]Next Steps:[/info] Phase 2 (ProteinMPNN Sequence Design) will be run next.")

if __name__ == "__main__":
    main()
