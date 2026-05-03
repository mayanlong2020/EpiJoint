import os, glob, pickle, shutil
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
from utils.logger import console
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

def get_plddt(pdb_file):
    with open(pdb_file, 'r') as f:
        bfactors = [float(line[60:66]) for line in f if line.startswith("ATOM  ") and line[12:16].strip() == "CA"]
    return np.mean(bfactors) if bfactors else 0.0

def run_phase4(run_dir: str):
    console.rule("[bold cyan]Phase 4: Structural QA & Unified Delivery (Enhanced)[/bold cyan]")
    
    af2_dir = os.path.join(run_dir, "03_validation", "03_af2_pdbs")
    mpnn_dir = os.path.join(run_dir, "02_mpnn_seqs", "seqs")
    backbone_dir = os.path.join(run_dir, "01_backbones")
    
    # 🎯 建立统一的交付目录
    delivery_dir = os.path.join(run_dir, "04_delivery")
    best_assets_dir = os.path.join(delivery_dir, "best_candidates")
    winner_dir = os.path.join(delivery_dir, "winner")
    os.makedirs(best_assets_dir, exist_ok=True)
    os.makedirs(winner_dir, exist_ok=True)
    
    native_pdb = os.path.abspath(os.path.join(run_dir, "../../inputs/rsvf_3_motifs.pdb"))
    parser = PDBParser(QUIET=True)
    native_struct = parser.get_structure("native", native_pdb)
    native_model = native_struct[0]
    
    # 预计算原始 Motif 集群中心（用于朝向判定）
    native_ca_atoms = [atom for atom in native_struct.get_atoms() if atom.get_name() == 'CA']
    native_center_all = np.mean([a.get_coord() for a in native_ca_atoms], axis=0)
    
    # 获取所有 Motif 链 ID
    native_chains = [chain.id for chain in native_model]
    console.print(f"[info]Detected {len(native_chains)} motifs in native structure: {native_chains}[/info]")
    
    results = []
    af2_pdbs = glob.glob(os.path.join(af2_dir, "*_unrelaxed_rank_001*.pdb"))
    console.print(f"[info]Evaluating exact Motif fidelity and orientation for {len(af2_pdbs)} structures...[/info]")
    
    for pdb in af2_pdbs:
        filename = os.path.basename(pdb)
        name = filename.split("_unrelaxed")[0] 
        base_design = name.split("_seq")[0]  
        
        trb_path = os.path.join(backbone_dir, f"{base_design}.trb")
        if not os.path.exists(trb_path): continue
            
        with open(trb_path, 'rb') as f: trb = pickle.load(f)
        plddt = get_plddt(pdb)
        
        try:
            pred_struct = parser.get_structure("pred", pdb)
            pred_model = pred_struct[0]
            pred_ca_all = [a for a in pred_model.get_atoms() if a.get_name() == 'CA']
            scaffold_center = np.mean([a.get_coord() for a in pred_ca_all], axis=0)
            
            res_data = {"Design_ID": name, "AF2_pLDDT": plddt}
            
            # 分 Motif 计算指标
            for chain_id in native_chains:
                native_motif_ca = [a for a in native_model[chain_id].get_atoms() if a.get_name() == 'CA']
                native_motif_center = np.mean([a.get_coord() for a in native_motif_ca], axis=0)
                # 原始向外向量
                native_out_vec = native_motif_center - native_center_all
                
                # 寻找预测结构中对应的原子
                pred_motif_ca = []
                for (ref_c, ref_r), (hal_c, hal_r) in zip(trb['con_ref_pdb_idx'], trb['con_hal_pdb_idx']):
                    if ref_c == chain_id:
                        try:
                            atom = pred_model[hal_c][(' ', int(hal_r), ' ')]['CA']
                            pred_motif_ca.append(atom)
                        except: continue
                
                if not pred_motif_ca or len(pred_motif_ca) != len(native_motif_ca): continue
                
                # 1. 内部 RMSD (只针对当前 motif)
                sup = Superimposer()
                sup.set_atoms(native_motif_ca, pred_motif_ca)
                motif_rmsd = sup.rms
                
                # 2. 朝向评分 (Orientation)
                rot, tran = sup.rotran
                pred_out_vec = np.dot(native_out_vec, rot) # 将原始向外向量旋转到当前坐标系
                
                motif_center_pred = np.mean([a.get_coord() for a in pred_motif_ca], axis=0)
                vec_to_motif = motif_center_pred - scaffold_center # 脚手架中心到 motif 的向量
                
                # 余弦相似度：1.0 表示完美朝外
                dot = np.dot(pred_out_vec, vec_to_motif) / (np.linalg.norm(pred_out_vec) * np.linalg.norm(vec_to_motif))
                orientation = np.clip(dot, -1.0, 1.0)
                
                # 3. 位阻评分 (Clash Check)
                clashes = 0
                for a in pred_ca_all:
                    if a in pred_motif_ca: continue
                    v = a.get_coord() - motif_center_pred
                    dist = np.linalg.norm(v)
                    if dist < 10.0: # 只检测 motif 周围 10A
                        # 如果原子位于“向外”向量的一侧，即视为可能遮挡结合面
                        if np.dot(v, pred_out_vec) > 0:
                            if dist < 4.5: clashes += 1 # 严格碰撞阈值
                
                res_data[f"Motif_{chain_id}_RMSD"] = motif_rmsd
                res_data[f"Motif_{chain_id}_Ori"] = orientation
                res_data[f"Motif_{chain_id}_Clash"] = clashes
            
            # 综合统计
            all_rmsds = [res_data[f"Motif_{c}_RMSD"] for c in native_chains if f"Motif_{c}_RMSD" in res_data]
            all_oris = [res_data[f"Motif_{c}_Ori"] for c in native_chains if f"Motif_{c}_Ori" in res_data]
            if not all_rmsds: continue
            
            res_data["Max_Motif_RMSD"] = max(all_rmsds)
            res_data["Min_Motif_Ori"] = min(all_oris) if all_oris else -1.0
            res_data["Total_Clashes"] = sum([res_data[f"Motif_{c}_Clash"] for c in native_chains if f"Motif_{c}_Clash" in res_data])
            
            res_data["PDB_Path"] = pdb
            res_data["TRB_Path"] = trb_path
            res_data["FASTA_Path"] = os.path.join(run_dir, "03_validation", "02_filtered_fastas_for_af2", f"{name}.fasta")
            
            results.append(res_data)
        except Exception as e:
            continue
            
    df = pd.DataFrame(results)
    if df.empty:
        console.print("[error]QA Failed. No valid structures could be aligned.[/error]")
        return
        
    # 排名逻辑：Max RMSD 越小越好，pLDDT 越高越好
    df_sorted = df.sort_values(by=['Max_Motif_RMSD', 'AF2_pLDDT'], ascending=[True, False])
    
    # 1. 导出全体质检报告
    report_path = os.path.join(delivery_dir, "QA_Report.csv")
    cols = ['Design_ID', 'AF2_pLDDT', 'Max_Motif_RMSD', 'Min_Motif_Ori', 'Total_Clashes'] + \
           [f"Motif_{c}_RMSD" for c in native_chains] + \
           [f"Motif_{c}_Ori" for c in native_chains]
    df_sorted[cols].to_csv(report_path, index=False)
    
    # 2. 动态筛选标准
    # 每个 Motif RMSD < 1.5 且每个 Motif Ori > 0 (朝向外侧)
    qualified = df_sorted[
        (df_sorted['AF2_pLDDT'] > 75) & 
        (df_sorted['Max_Motif_RMSD'] < 1.5) &
        (df_sorted['Min_Motif_Ori'] > 0.0) & 
        (df_sorted['Total_Clashes'] < 5)
    ]
    
    # 清理旧的交付资产，确保一致性
    if os.path.exists(best_assets_dir): shutil.rmtree(best_assets_dir)
    if os.path.exists(winner_dir): shutil.rmtree(winner_dir)
    os.makedirs(best_assets_dir, exist_ok=True)
    os.makedirs(winner_dir, exist_ok=True)
    
    # 如果没有设计完全达标，则放宽标准取前 10 名放入 best_candidates，但 winner 保持为空
    if qualified.empty:
        console.print("[warning]No candidates met strict criteria (All Motifs RMSD < 1.5 & Ori > 0). Listing top 10 in best_candidates only.[/warning]")
        display_candidates = df_sorted.head(10)
    else:
        display_candidates = qualified
        console.print(f"\n[bold gold1]✨ Found {len(qualified)} WINNERS meeting all criteria! ✨[/bold gold1]")
    
    for i, (idx, row) in enumerate(display_candidates.iterrows()):
        design_id = row['Design_ID']
        
        # 复制到合格目录
        shutil.copy(row['PDB_Path'], os.path.join(best_assets_dir, f"{design_id}.pdb"))
        shutil.copy(row['FASTA_Path'], os.path.join(best_assets_dir, f"{design_id}.fasta"))
        shutil.copy(row['TRB_Path'], os.path.join(best_assets_dir, f"{design_id}.trb"))
        
        # 如果是合格的 Winner，则进入 winner 目录
        if not qualified.empty and design_id in qualified['Design_ID'].values:
            shutil.copy(row['PDB_Path'], os.path.join(winner_dir, f"{design_id}.pdb"))
            shutil.copy(row['FASTA_Path'], os.path.join(winner_dir, f"{design_id}.fasta"))
            console.print(f"  🏆 [bold gold1]WINNER[/bold gold1]: [highlight]{design_id}[/highlight] | Max RMSD: {row['Max_Motif_RMSD']:.2f} | Min Ori: {row['Min_Motif_Ori']:.2f}")
        else:
            console.print(f"  ► {design_id} | Max RMSD: {row['Max_Motif_RMSD']:.2f} | Min Ori: {row['Min_Motif_Ori']:.2f} | pLDDT: {row['AF2_pLDDT']:.1f}")

    console.print(f"\n[success]✔ Pipeline Complete! {len(display_candidates)} delivery assets generated.[/success]")
    console.print(f"📦 Check your payload at: {delivery_dir}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        run_phase4(sys.argv[1])
