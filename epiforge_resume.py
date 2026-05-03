import os
import time
from utils.logger import console

# 导入所有独立的流水线模块
from scripts.phase0_prepare import run_phase0
from scripts.phase1_generate import run_phase1
from scripts.phase2_mpnn import run_phase2
from scripts.phase3_validate import run_phase3
from scripts.phase4_qa import run_phase4

# ==========================================
# ⚙️ 恢复运行配置 (应与 epiforge_main.py 保持一致)
# ==========================================
NUM_DESIGNS = 1         
CONCURRENCY_CPU = 14       
# ==========================================

def main():
    console.print("\n[bold highlight]=================================================[/bold highlight]")
    console.print("[bold highlight]   🧬 EpiForge Pipeline: Resuming Session...   [/bold highlight]")
    console.print("[bold highlight]=================================================[/bold highlight]\n")
    
    base_dir = os.path.abspath(os.path.dirname(__file__))
    inputs_dir = os.path.join(base_dir, "inputs")
    runs_dir = os.path.join(base_dir, "runs")
    
    # 自动探测最新的运行目录
    all_runs = sorted([os.path.join(runs_dir, d) for d in os.listdir(runs_dir) if d.startswith("run_")])
    if not all_runs:
        console.print("[error]No run directories found![/error]")
        return
        
    current_run_dir = all_runs[-1]
    console.print(f"[info]Latest Session Detected:[/info] {current_run_dir}")
    console.print(f"[info]Target Designs:[/info] {NUM_DESIGNS}\n")
    
    pipeline_start_time = time.time()
    
    try:
        # 🔴 Phase 0: 准备天然靶标 (通常是跳过)
        motif_pdb = run_phase0(inputs_dir)
        
        # 🟠 Phase 1: 拓扑骨架生成 (RFjoint2) - 现在支持断点续传
        run_phase1(current_run_dir, motif_pdb, num_designs=NUM_DESIGNS, max_workers=CONCURRENCY_CPU)
        
        # 🟡 Phase 2: 序列设计 (ProteinMPNN) - 现在支持断点续传
        run_phase2(current_run_dir)
        
        # 🟢 Phase 3: 结构验证 (ESMFold + AF2) - 现在支持断点续传
        run_phase3(current_run_dir)
        
        # 🔵 Phase 4: 交付打包
        run_phase4(current_run_dir)
        
    except Exception as e:
        console.print_exception(show_locals=True)
        console.print(f"\n[bold red]✖ Pipeline resume failed.[/bold red]")
        return

    total_time = time.time() - pipeline_start_time
    console.rule("[bold highlight]🎉 Pipeline Resumed and Completed! 🎉[/bold highlight]")
    console.print(f"📦 [bold]All assets perfectly packaged at:[/bold] {current_run_dir}/04_delivery")

if __name__ == "__main__":
    main()
