import os
import time
from datetime import datetime
from utils.logger import console

# 导入所有独立的流水线模块
from scripts.phase0_prepare import run_phase0
from scripts.phase1_generate import run_phase1
from scripts.phase2_mpnn import run_phase2
from scripts.phase3_validate import run_phase3
from scripts.phase4_qa import run_phase4
from scripts.phase5_binding import run_phase5

# ==========================================
# ⚙️ 核心运行配置 (在这里修改你的目标策略)
# ==========================================
NUM_DESIGNS = 1         # 🚀 你想生成的骨架数量（今晚通宵建议改为 5000）
CONCURRENCY_CPU = 14       # M4 Max 的性能核+能效核总数（用于 Phase 1）
# ==========================================

def main():
    console.print("\n[bold highlight]=================================================[/bold highlight]")
    console.print("[bold highlight] 🚀 EpiForge Pipeline: Fully Automated End-to-End 🚀 [/bold highlight]")
    console.print("[bold highlight]=================================================[/bold highlight]\n")
    
    # 1. 基础架构初始化
    base_dir = os.path.abspath(os.path.dirname(__file__))
    inputs_dir = os.path.join(base_dir, "inputs")
    runs_dir = os.path.join(base_dir, "runs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(runs_dir, exist_ok=True)
    
    # 创建带时间戳的全新运行隔离目录
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    current_run_dir = os.path.join(runs_dir, f"run_{timestamp}")
    os.makedirs(current_run_dir, exist_ok=True)
    
    console.print(f"[info]Initialized Master Session:[/info] {current_run_dir}")
    console.print(f"[info]Target Designs:[/info] {NUM_DESIGNS}\n")
    
    pipeline_start_time = time.time()
    
    try:
        # -----------------------------------------------------
        # 🔴 Phase 0: 准备天然靶标
        # -----------------------------------------------------
        motif_pdb = run_phase0(inputs_dir)
        
        # -----------------------------------------------------
        # 🟠 Phase 1: 拓扑骨架生成 (RFjoint2)
        # -----------------------------------------------------
        run_phase1(current_run_dir, motif_pdb, num_designs=NUM_DESIGNS, max_workers=CONCURRENCY_CPU)
        
        # -----------------------------------------------------
        # 🟡 Phase 2: 序列设计与约束 (ProteinMPNN)
        # -----------------------------------------------------
        run_phase2(current_run_dir)
        
        # -----------------------------------------------------
        # 🟢 Phase 3: 双重结构验证 (ESMFold 快筛 + AlphaFold2 精筛)
        # -----------------------------------------------------
        run_phase3(current_run_dir)
        
        # -----------------------------------------------------
        # 🔵 Phase 4: 3D 几何对齐 QA 与交付打包
        # -----------------------------------------------------
        run_phase4(current_run_dir)
        
        # -----------------------------------------------------
        # 🟣 Phase 5: 抗体抗原结合面测试 (ColabFold Multimer)
        # -----------------------------------------------------
        run_phase5(current_run_dir)
        
    except Exception as e:
        console.print_exception(show_locals=True)
        console.print(f"\n[bold red]✖ Pipeline crashed due to an unexpected error.[/bold red]")
        return

    # 汇总计算时间
    total_time = time.time() - pipeline_start_time
    hours, rem = divmod(total_time, 3600)
    minutes, seconds = divmod(rem, 60)
    
    console.rule("[bold highlight]🎉 EpiForge Pipeline Fully Completed! 🎉[/bold highlight]")
    console.print(f"[success]Total Execution Time:[/success] {int(hours)}h {int(minutes)}m {seconds:.1f}s")
    console.print(f"📦 [bold]All assets perfectly packaged at:[/bold] {current_run_dir}/04_delivery")

if __name__ == "__main__":
    main()
