import argparse
import os
import torch
import warnings

# 固化环境变量
os.environ["PYTORCH_ENABLE_MPS_FALLBACK"] = "1"
warnings.filterwarnings("ignore")

from transformers import AutoTokenizer, EsmForProteinFolding

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", type=str, required=True)
    parser.add_argument("--out_dir", type=str, required=True)
    parser.add_argument("--out_csv", type=str, required=True)
    parser.add_argument("--chunk_id", type=str, default="0")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    
    tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
    model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
    
    if torch.backends.mps.is_available():
        device = torch.device("mps")
    else:
        device = torch.device("cpu")
        
    model = model.to(device)
    model.eval()

    with open(args.fasta, "r") as f:
        lines = f.read().splitlines()
        
    names, seqs =[],[]
    for line in lines:
        if line.startswith(">"): names.append(line[1:])
        elif line.strip():
            if len(names) > len(seqs): seqs.append(line.strip())
            else: seqs[-1] += line.strip()

    results =[]
    with torch.no_grad():
        for i, (name, seq) in enumerate(zip(names, seqs)):
            out_pdb = os.path.join(args.out_dir, f"{name}.pdb")
            
            # 断点续传 & 修复 0-1 到 0-100 的数值映射
            if os.path.exists(out_pdb):
                with open(out_pdb, 'r') as f:
                    b_factors = [float(line[60:66]) for line in f if line.startswith("ATOM  ") and line[12:16].strip() == "CA"]
                if b_factors:
                    plddt = sum(b_factors) / len(b_factors)
                    if plddt <= 1.0: plddt *= 100  # 修复由于 HuggingFace 范围导致的偏低问题
                    results.append((name, plddt))
                    print(f"[Chunk {args.chunk_id}] {name} | pLDDT: {plddt:.2f} [Cached]")
                    continue
            
            # 推理
            inputs = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
            inputs = {k: v.to(device) for k, v in inputs.items()}
            outputs = model(**inputs)
            
            plddt = outputs.plddt[0, :len(seq)].mean().item() * 100
            pdb_content = model.output_to_pdb(outputs)[0]
            
            temp_pdb = out_pdb + f".tmp_{args.chunk_id}"
            with open(temp_pdb, "w") as f: f.write(pdb_content)
            os.rename(temp_pdb, out_pdb)
            
            results.append((name, plddt))
            print(f"[Chunk {args.chunk_id}] {name} | pLDDT: {plddt:.2f} [Computed]")
            
            # 主动清理显存碎片
            if str(device) == "mps":
                torch.mps.empty_cache()

    with open(args.out_csv, 'w') as f:
        f.write("design_name,esm_plddt\n")
        for name, plddt in results:
            f.write(f"{name},{plddt:.2f}\n")

if __name__ == "__main__":
    main()
