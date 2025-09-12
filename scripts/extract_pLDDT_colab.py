import os
import re
import sys
import pandas as pd

def ask(msg, default=None, allow_empty=False):
    tip = f" [{default}]" if default else ""
    while True:
        s = input(f"{msg}{tip}: ").strip()
        if not s and default: return default
        if s or allow_empty: return s
        print("Please enter a value.")

def ask_fasta_len_or_path(msg, default=""):
    s = input(f"{msg} [{default}]: ").strip()
    if not s and default: s = default
    if not s: return None
    if s.isdigit(): return int(s)
    if not os.path.isfile(s):
        print(f"FASTA not found: {s}")
        return None
    length = 0
    with open(s, "r") as f:
        for line in f:
            if line.startswith(">"):
                if length > 0: break
                else: continue
            length += len(line.strip())
    return length if length>0 else None

def parse_loop_ranges(s):
    if not s.strip(): return []
    out=[]
    for part in [p.strip() for p in s.split(";") if p.strip()]:
        m = re.match(r"^\s*(-?\d+)\s*[-:]\s*(-?\d+)\s*$", part)
        if m:
            a,b=int(m.group(1)),int(m.group(2))
            if a>b: a,b=b,a
            out.append((a,b))
    return out

def detect_state(name):
    nl=name.lower()
    if "unrelaxed" in nl: return "unrelaxed"
    if "relaxed" in nl: return "relaxed"
    return "unknown"

def find_pdbs(root):
    out=[]
    for base,_,files in os.walk(root):
        for fn in files:
            if fn.lower().endswith(".pdb"): out.append(os.path.join(base,fn))
    return sorted(out)

def parse_ca_plddt(pdb):
    vals=[]; resmap={}
    with open(pdb) as f:
        for line in f:
            if not line.startswith("ATOM"): continue
            if line[12:16].strip()!="CA": continue
            res=line[22:26].strip()
            try: res=int(res)
            except: continue
            try: val=float(line[60:66].strip())
            except: continue
            vals.append(val); resmap[res]=val
    return vals,resmap

def summarize(vals,thr70=70,thr90=90):
    if not vals: return 0,None,None,None
    n=len(vals); mean=round(sum(vals)/n,2)
    pct70=round(sum(v>=thr70 for v in vals)/n*100,1)
    pct90=round(sum(v>=thr90 for v in vals)/n*100,1)
    return n,mean,pct70,pct90

def loop_stats(resmap,ranges):
    out={}
    for a,b in ranges:
        vals=[resmap[r] for r in range(a,b+1) if r in resmap]
        tag=f"Loop_{a}-{b}"
        if vals:
            out[f"{tag}_mean"]=round(sum(vals)/len(vals),2)
            out[f"{tag}_min"]=round(min(vals),2)
            out[f"{tag}_max"]=round(max(vals),2)
            out[f"{tag}_n"]=len(vals)
        else:
            out[f"{tag}_mean"]=None; out[f"{tag}_min"]=None; out[f"{tag}_max"]=None; out[f"{tag}_n"]=0
    return out

def pick_best(df):
    df["__relaxed__"]=(df["State"]=="relaxed").astype(int)
    best_idx=[]
    for enz,sub in df.groupby("Enzyme",sort=False):
        sub2=sub.sort_values(by=["__relaxed__","Mean_pLDDT"],ascending=[False,False])
        best_idx.append(sub2.index[0])
    df["★"]=""
    df.loc[best_idx,"★"]="★"
    df.drop(columns=["__relaxed__"],inplace=True)
    return df

# -------- run --------

root=ask("ColabFold output directory (absolute path)","/path/to/colabfold_outputs")
if not os.path.isdir(root): sys.exit(f"Directory not found: {root}")
name=ask("Name (used for Enzyme & CSV prefix)","BfIMTD")
fastalen=ask_fasta_len_or_path("FASTA length (int) OR FASTA file path (Enter to skip)","")
loopranges=parse_loop_ranges(ask("Loop ranges e.g. 392-402;450-460 (Enter to skip)","",allow_empty=True))

pdbs=find_pdbs(root)
if not pdbs: sys.exit("No PDBs found.")

records=[]
for pdb in pdbs:
    vals,resmap=parse_ca_plddt(pdb)
    n,mean,p70,p90=summarize(vals)
    row={"Enzyme":name,"ModelFileName":os.path.basename(pdb),
         "Mean_pLDDT":mean,"Residues_gt70_%":p70,"Residues_gt90_%":p90,
         "Length_CA":n,"FASTA_length":fastalen,"State":detect_state(pdb)}
    if loopranges: row.update(loop_stats(resmap,loopranges))
    records.append(row)

df=pd.DataFrame(records)
df=pick_best(df)

# save to sibling "summary"
save_dir=os.path.join(os.path.dirname(root),"summary")
os.makedirs(save_dir,exist_ok=True)
out=os.path.join(save_dir,f"{name}_plddt_summary.csv")
df.to_csv(out,index=False)

print(f"✅ Saved CSV: {out}")