import os
import sys
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import torch
    import esm
except ImportError:
    print("ERROR_FAIR_ESM_MISSING")
    sys.exit(1)

# Worker 进程全局变量
worker_model = None

def init_worker():
    global worker_model
    # 限制单进程的线程数，防止 4 个进程互相踩踏 CPU 核心
    torch.set_num_threads(3)
    os.environ["OMP_NUM_THREADS"] = "3"
    os.environ["MKL_NUM_THREADS"] = "3"
    
    # 将模型加载到 CPU
    worker_model = esm.pretrained.esmfold_v1().eval().cpu()

def fold_sequence(task_args):
    fasta_path, out_pdb_dir = task_args
    name = os.path.basename(fasta_path).replace('.fasta', '')
    out_pdb = os.path.join(out_pdb_dir, f"{name}.pdb")

    # 🚀 避免重复计算：如果 PDB 已存在，直接解析现有的 pLDDT（利用 B-factor 列）
    if os.path.exists(out_pdb):
        with open(out_pdb, 'r') as f:
            b_factors = [float(line[60:66]) for line in f if line.startswith("ATOM  ") and line[12:16].strip() == "CA"]
        if b_factors:
            return name, sum(b_factors)/len(b_factors), True # True 表示复用缓存

    with open(fasta_path, 'r') as f:
        lines = f.read().splitlines()
        seq = lines[1] if len(lines) > 1 else ""

    if not seq: return name, 0.0, False

    with torch.no_grad():
        output = worker_model.infer(seq)
        pdb_str = worker_model.output_to_pdb(output)[0]
        plddt = output['mean_plddt'][0].item() * 100

    # 🚀 避免 IO 竞态竞争：先写临时文件，然后原子化重命名
    temp_pdb = out_pdb + f".tmp_{os.getpid()}"
    with open(temp_pdb, 'w') as f:
        f.write(pdb_str)
    os.rename(temp_pdb, out_pdb)

    return name, plddt, False

def main():
    fasta_dir = sys.argv[1]
    out_pdb_dir = sys.argv[2]
    out_csv = sys.argv[3]
    concurrency = int(sys.argv[4]) if len(sys.argv) > 4 else 4

    os.makedirs(out_pdb_dir, exist_ok=True)
    fastas = glob.glob(os.path.join(fasta_dir, "*.fasta"))

    print("[ESMFold] Pre-loading weights safely to avoid cache race conditions...")
    _ = esm.pretrained.esmfold_v1() # 预触发下载验证
    del _

    print(f"\n[ESMFold] Igniting parallel inference with {concurrency} workers on {len(fastas)} sequences...")

    results =[]
    tasks = [(f, out_pdb_dir) for f in fastas]

    # 启动多进程池
    with ProcessPoolExecutor(max_workers=concurrency, initializer=init_worker) as executor:
        futures =[executor.submit(fold_sequence, t) for t in tasks]
        for i, future in enumerate(as_completed(futures)):
            name, plddt, cached = future.result()
            results.append((name, plddt))
            status = "[Cached]" if cached else "[Computed]"
            print(f"[ESMFold] {i+1}/{len(fastas)} | {name} | pLDDT: {plddt:.2f} {status}", flush=True)

    # 汇总写入 CSV
    with open(out_csv, 'w') as f:
        f.write("design_name,esm_plddt\n")
        for name, plddt in sorted(results, key=lambda x: x[1], reverse=True):
            f.write(f"{name},{plddt:.2f}\n")

    print(f"\n[ESMFold] Complete! Saved to {out_csv}")

if __name__ == "__main__":
    main()
