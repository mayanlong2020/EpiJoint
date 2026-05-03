import os
import glob
import json
import pickle
import subprocess
import shutil
import time
from utils.logger import console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn, TimeElapsedColumn

# ProteinMPNN 路径配置
MPNN_DIR = os.path.expanduser("~/Development/BioDesign/ProteinMPNN")
HELPER_DIR = os.path.join(MPNN_DIR, "helper_scripts")

def run_phase2(run_dir: str):
    console.rule("[bold cyan]Phase 2: Sequence Design (ProteinMPNN - High Performance Batch)[/bold cyan]")
    backbone_dir = os.path.join(run_dir, "01_backbones")
    out_dir = os.path.join(run_dir, "02_mpnn_seqs")
    temp_pdb_dir = os.path.join(out_dir, "temp_pdbs_to_process")
    seqs_out_dir = os.path.join(out_dir, "seqs")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(seqs_out_dir, exist_ok=True)
    
    # 1. 筛选需要处理的设计 (Resume 逻辑)
    all_pdbs = glob.glob(os.path.join(backbone_dir, "*_0.pdb"))
    already_done_count = len(glob.glob(os.path.join(seqs_out_dir, "*.fa")))
    
    pdbs_to_process = []
    all_fixed_positions = {}

    for pdb_path in all_pdbs:
        base_name = os.path.basename(pdb_path).replace(".pdb", "")
        mpnn_out_seqs = os.path.join(seqs_out_dir, f"{base_name}.fa")
        
        if os.path.exists(mpnn_out_seqs):
            continue # 已完成，跳过
            
        trb_file = pdb_path.replace(".pdb", ".trb")
        if not os.path.exists(trb_file): continue
        with open(trb_file, 'rb') as f: trb = pickle.load(f)
        
        fixed_pos = [int(idx) for chain, idx in trb.get('con_hal_pdb_idx', []) if chain == 'A']
        if fixed_pos:
            pdbs_to_process.append(pdb_path)
            all_fixed_positions[base_name] = {"A": fixed_pos}

    if not pdbs_to_process:
        console.print(f"[success]✔ All 200 sequences already designed. Skipping Phase 2.[/success]")
        return

    console.print(f"[info]Already finished: {already_done_count} | Pending: {len(pdbs_to_process)}[/info]")
    console.print(f"[info]Preparing batch processing for the remaining {len(pdbs_to_process)} designs...[/info]")

    # 2. 预处理
    if os.path.exists(temp_pdb_dir): shutil.rmtree(temp_pdb_dir)
    os.makedirs(temp_pdb_dir)
    
    unified_jsonl = os.path.join(out_dir, "combined_fixed_positions.jsonl")
    parsed_jsonl = os.path.join(out_dir, "parsed_pdbs.jsonl")
    
    with open(unified_jsonl, 'w') as f:
        json.dump(all_fixed_positions, f)
        f.write("\n")

    for name in all_fixed_positions.keys():
        os.symlink(os.path.join(backbone_dir, f"{name}.pdb"), os.path.join(temp_pdb_dir, f"{name}.pdb"))

    # 第一步：解析坐标
    parse_cmd = [
        "conda", "run", "-n", "proteinmpnn",
        "python", os.path.join(HELPER_DIR, "parse_multiple_chains.py"),
        "--input_path", temp_pdb_dir,
        "--output_path", parsed_jsonl
    ]
    
    try:
        with console.status("[bold info]Pre-parsing PDB coordinates...[/bold info]"):
            subprocess.run(parse_cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        console.print(f"[error]PDB Parsing failed![/error]")
        return

    # 第二步：批量设计 + 进度监控
    design_cmd = [
        "conda", "run", "-n", "proteinmpnn",
        "python", os.path.join(MPNN_DIR, "protein_mpnn_run.py"),
        "--jsonl_path", parsed_jsonl,
        "--fixed_positions_jsonl", unified_jsonl,
        "--out_folder", out_dir,
        "--num_seq_per_target", "16",
        "--sampling_temp", "0.1",
        "--batch_size", "8",
        "--suppress_print", "1"
    ]
    
    total_to_do = len(pdbs_to_process)
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        BarColumn(),
        TaskProgressColumn(),
        TimeElapsedColumn(),
        console=console
    ) as progress:
        task = progress.add_task(f"[cyan]Designing {total_to_do} sequences in bulk...", total=total_to_do)
        
        # 异步启动
        process = subprocess.Popen(design_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        
        while process.poll() is None:
            # 实时扫描输出目录，计算增长量
            current_count = len(glob.glob(os.path.join(seqs_out_dir, "*.fa")))
            newly_done = current_count - already_done_count
            progress.update(task, completed=min(newly_done, total_to_do))
            time.sleep(1) # 每秒刷新一次
            
        # 最终确认状态
        if process.returncode != 0:
            console.print(f"[error]ProteinMPNN Design failed with exit code {process.returncode}[/error]")
        else:
            progress.update(task, completed=total_to_do)

    # 4. 清理
    if os.path.exists(temp_pdb_dir): shutil.rmtree(temp_pdb_dir)
    if os.path.exists(parsed_jsonl): os.remove(parsed_jsonl)

    final_seq_count = len(glob.glob(os.path.join(seqs_out_dir, "*.fa")))
    console.print(f"\n[success]✔ Phase 2 Completed! Current total: {final_seq_count}/200 designs finished.[/success]")
