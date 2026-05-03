# VER: BOLTZ_V10 (Conda-Native Patch)
import os, glob, shutil, subprocess, json, csv, sys, yaml
from utils.logger import console

def extract_fasta_records(fasta_path):
    records = []
    curr_id = ""
    curr_seq = ""
    if not os.path.exists(fasta_path): return []
    with open(fasta_path, 'r') as f:
        for line in f:
            if line.startswith(">"):
                if curr_id and curr_seq:
                    records.append({"id": curr_id, "seq": curr_seq})
                curr_id = line[1:].strip().split()[0]
                curr_seq = ""
            else:
                curr_seq += line.strip()
    if curr_id and curr_seq:
        records.append({"id": curr_id, "seq": curr_seq})
    return records

def run_phase5(run_dir: str):
    console.rule("[bold cyan]Phase 5: Antibody-Antigen Binding Analysis (Boltz-2 Conda-Native V10)[/bold cyan]")
    
    delivery_dir = os.path.join(run_dir, "04_delivery")
    winner_dir = os.path.join(delivery_dir, "winner")
    binding_dir = os.path.join(run_dir, "05_binding_analysis")
    yaml_prep_dir = os.path.join(binding_dir, "boltz_inputs")
    results_dir = os.path.join(binding_dir, "boltz_results")
    patch_dir = os.path.join(binding_dir, "patch")
    
    os.makedirs(yaml_prep_dir, exist_ok=True)
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(patch_dir, exist_ok=True)
    
    # sitecustomize.py 猴子补丁
    patch_path = os.path.join(patch_dir, "sitecustomize.py")
    with open(patch_path, "w") as f:
        f.write("import torch\n")
        f.write("import functools\n")
        f.write("_orig_load = torch.load\n")
        f.write("@functools.wraps(_orig_load)\n")
        f.write("def _patched_load(*args, **kwargs):\n")
        f.write("    if 'weights_only' in kwargs: kwargs['weights_only'] = False\n")
        f.write("    return _orig_load(*args, **kwargs)\n")
        f.write("torch.load = _patched_load\n")
        f.write("print('  [Conda-Patch] Torch.load bypass active.')\n")
    
    boltz_cache = os.path.expanduser("~/Development/BioDesign/Boltz_2/boltz_cache")
    
    # 1. 识别获胜者
    winner_fastas = glob.glob(os.path.join(winner_dir, "design_*.fasta"))
    if not winner_fastas:
        alt = os.path.join(winner_dir, "Target_Antigen_Winner.fasta")
        if os.path.exists(alt): winner_fastas = [alt]
        else:
            console.print("[warning]No winner fastas found.[/warning]")
            return

    # 2. 识别抗体
    antibodies = {
        "SiteII": os.path.join(winner_dir, "Antibody_Site_II_Palivizumab.fasta"),
        "SiteIV": os.path.join(winner_dir, "Antibody_Site_IV_101F.fasta"),
        "SiteV":  os.path.join(winner_dir, "Antibody_Site_V_hRSV90.fasta")
    }
    
    # 3. 生成 YAML
    tasks = []
    for winner_fa in winner_fastas:
        w_records = extract_fasta_records(winner_fa)
        if not w_records: continue
        w_name = os.path.basename(winner_fa).replace(".fasta", "")
        for ab_site, ab_fa in antibodies.items():
            ab_records = extract_fasta_records(ab_fa)
            if not ab_records: continue
            task_name = f"{w_name}_vs_{ab_site}"
            yaml_path = os.path.join(yaml_prep_dir, f"{task_name}.yaml")
            boltz_data = {"version": 1, "sequences": []}
            boltz_data["sequences"].append({"protein": {"id": "A", "sequence": w_records[0]["seq"], "msa": "empty"}})
            for i, rec in enumerate(ab_records):
                chain_id = "H" if i == 0 else "L"
                boltz_data["sequences"].append({"protein": {"id": chain_id, "sequence": rec["seq"], "msa": "empty"}})
            with open(yaml_path, 'w') as f: yaml.dump(boltz_data, f)
            tasks.append(task_name)

    # 4. 执行
    env = os.environ.copy()
    env["BOLTZ_CACHE_DIR"] = boltz_cache
    # 强制注入补丁路径
    env["PYTHONPATH"] = f"{patch_dir}:{env.get('PYTHONPATH', '')}"
    env["KMP_DUPLICATE_LIB_OK"] = "TRUE"
    
    for task_name in tasks:
        yaml_path = os.path.join(yaml_prep_dir, f"{task_name}.yaml")
        task_out_dir = os.path.join(results_dir, task_name)
        
        pred_base = os.path.join(task_out_dir, f"boltz_results_{task_name}", "predictions", task_name)
        if glob.glob(os.path.join(pred_base, "confidence_*.json")):
            console.print(f"  ⏭  Already done: {task_name}")
            continue

        console.print(f"  🚀 [bold]Boltz-2 Predicting (Conda Run):[/bold] {task_name}")
        # 使用 conda run 确保完全进入 boltz2 环境
        cmd = [
            "conda", "run", "--no-capture-output", "-n", "boltz2",
            "boltz", "predict", yaml_path,
            "--out_dir", task_out_dir,
            "--cache", boltz_cache,
            "--accelerator", "gpu",
            "--devices", "1",
            "--output_format", "pdb"
        ]
        try:
            subprocess.run(cmd, env=env, check=True)
            console.print(f"  ✅ [bold green]Success:[/bold green] {task_name}")
        except Exception as e:
            console.print(f"  ❌ Boltz-2 failed: {e}")
            continue

    # 5. 汇总
    console.print("\n[info]Generating Final Binding Report...[/info]")
    report_rows = []
    for task in tasks:
        pred_base = os.path.join(results_dir, task, f"boltz_results_{task}", "predictions", task)
        conf_files = glob.glob(os.path.join(pred_base, "confidence_*.json"))
        if conf_files:
            try:
                with open(conf_files[0], 'r') as f:
                    data = json.load(f)
                report_rows.append({
                    "Design": task.split("_vs_")[0],
                    "Antibody": task.split("_vs_")[1],
                    "Score": round(data.get("confidence_score", 0.0), 3),
                    "ipTM": round(data.get("iptm", 0.0), 3)
                })
                pdb_files = glob.glob(os.path.join(pred_base, "*.pdb"))
                if pdb_files:
                    shutil.copy(pdb_files[0], os.path.join(delivery_dir, f"Binding_{task}.pdb"))
            except Exception as e:
                console.print(f"[error]Parse error {task}: {e}[/error]")

    if report_rows:
        report_path = os.path.join(delivery_dir, "Binding_Analysis_Report.md")
        with open(report_path, "w") as f:
            f.write("# 抗体抗原结合分析报告 (Boltz-2 V10)\n\n| 设计方案 | 测试抗体位点 | Confidence Score | ipTM |\n| :--- | :--- | :--- | :--- |\n")
            for row in report_rows:
                f.write(f"| {row['Design']} | {row['Antibody']} | {row['Score']} | {row['ipTM']} |\n")
        console.print(f"[success]✔ Final report generated at {report_path}[/success]")

if __name__ == "__main__":
    if len(sys.argv) > 1: run_phase5(sys.argv[1])
