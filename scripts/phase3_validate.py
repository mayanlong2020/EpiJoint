import os, glob, shutil, subprocess
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from utils.logger import console

# 固化环境变量
os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"

def chunk_dict(d, chunk_size=50):
    it = iter(d)
    for i in range(0, len(d), chunk_size):
        yield {k: d[k] for k in[next(it) for _ in range(min(chunk_size, len(d) - i))]}

def run_phase3(run_dir: str):
    console.rule("[bold cyan]Phase 3: High-Concurrency Dual-Stage Validation[/bold cyan]")
    mpnn_dir = os.path.join(run_dir, "02_mpnn_seqs", "seqs")
    out_dir = os.path.join(run_dir, "03_validation")
    
    esm_pdb_dir = os.path.join(out_dir, "01_esmfold_pdbs")
    chunks_dir = os.path.join(out_dir, "01_fasta_chunks")
    af2_fasta_dir = os.path.join(out_dir, "02_filtered_fastas_for_af2")
    af2_out_dir = os.path.join(out_dir, "03_af2_pdbs")
    
    for d in[esm_pdb_dir, chunks_dir, af2_fasta_dir, af2_out_dir]:
        os.makedirs(d, exist_ok=True)
    
    # 1. 提取序列并划分为 Chunks
    seqs = {}
    for fa in glob.glob(os.path.join(mpnn_dir, "*.fa")):
        base = os.path.basename(fa).replace(".fa", "")
        with open(fa, 'r') as f: lines = f.read().splitlines()
        curr = ""
        for line in lines:
            if line.startswith(">"):
                curr = f"{base}_seq{line.split('sample=')[1].split(',')[0]}" if "sample=" in line else "native"
            elif curr != "native" and curr:
                seqs[curr] = line
                
    total_seqs = len(seqs)
    console.print(f"[success]✔ Loaded {total_seqs} sequences.[/success]")
    
    chunk_files =[]
    for i, chunk in enumerate(chunk_dict(seqs, chunk_size=50)):
        chunk_fa = os.path.join(chunks_dir, f"chunk_{i}.fasta")
        with open(chunk_fa, "w") as f:
            for name, seq in chunk.items():
                f.write(f">{name}\n{seq}\n")
        chunk_files.append(chunk_fa)
        
    console.print(f"[info]Split into {len(chunk_files)} chunks (50 seqs/chunk) to prevent memory leaks.[/info]\n")
    
    # 2. 多并发调度 ESMFold Worker
    console.print("[info]Igniting Concurrent ESMFold GPU Workers...[/info]")
    def run_worker(chunk_fa, chunk_id):
        out_csv = os.path.join(chunks_dir, f"results_{chunk_id}.csv")
        # 💡 断点续传逻辑：如果该 Chunk 已有结果，直接返回
        if os.path.exists(out_csv):
            return out_csv
            
        cmd =[
            "conda", "run", "--no-capture-output", "-n", "esmfold_env",
            "python", "scripts/esmfold_hf.py",
            "--fasta", chunk_fa,
            "--out_dir", esm_pdb_dir,
            "--out_csv", out_csv,
            "--chunk_id", str(chunk_id)
        ]
        # 继承主进程的环境变量
        subprocess.run(cmd, env=os.environ.copy())
        return out_csv

    # 并发数为 4，M4 Max的 GPU 可以完美处理 4 路流
    csv_files =[]
    with ThreadPoolExecutor(max_workers=4) as executor:
        futures = {executor.submit(run_worker, f, i): i for i, f in enumerate(chunk_files)}
        for future in as_completed(futures):
            csv_files.append(future.result())
            
    # 3. 合并打分并过滤
    console.print("\n[info]Merging scores and applying pLDDT threshold...[/info]")
    dfs =[pd.read_csv(csv) for csv in csv_files if os.path.exists(csv)]
    if not dfs:
        console.print("[error]No CSV results generated! Halting.[/error]")
        return
        
    df = pd.concat(dfs).sort_values(by="esm_plddt", ascending=False)
    final_csv = os.path.join(out_dir, "esmfold_scores.csv")
    df.to_csv(final_csv, index=False)
    
    # 设定合理阈值，保留精华
    THRESHOLD = 60.0
    passed_df = df[df['esm_plddt'] >= THRESHOLD]
    passed_count = len(passed_df)
    console.print(f"[warning]Kept {passed_count}/{total_seqs} sequences with ESM pLDDT >= {THRESHOLD}[/warning]")
    
    if passed_count == 0:
        console.print("[error]No sequences passed. Halting.[/error]")
        return
        
    for name in passed_df['design_name']:
        with open(os.path.join(af2_fasta_dir, f"{name}.fasta"), "w") as f:
            f.write(f">{name}\n{seqs[name]}\n")
            
    # 4. AF2 精确渲染
    console.print(f"\n[info]Igniting AlphaFold2 High-Fidelity Validation on {passed_count} candidates...[/info]")
    af2_cmd =[
        "conda", "run", "--no-capture-output", "-n", "colabfold_m4",
        "colabfold_batch", af2_fasta_dir, af2_out_dir,
        "--data", os.path.expanduser("~/Development/BioDesign/ColabFold/models"), 
        "--msa-mode", "single_sequence", "--model-type", "alphafold2_ptm",
        "--num-models", "1", "--num-recycle", "3", "--disable-unified-memory"       
    ]
    subprocess.run(af2_cmd, env=os.environ.copy())
    console.print(f"\n[success]✔ Phase 3 Completed![/success]")
