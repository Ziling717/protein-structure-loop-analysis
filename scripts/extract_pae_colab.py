# ====== Minimal PAE analysis for ColabFold (save to sibling "summary/") ======
# What this does:
# - Find predicted_aligned_error*.json in the given ColabFold output folder (prefer rank_001/1)
# - Compute mean PAE (Å) for: Global, Loop-Loop, Loop-Core (sym.), Core-Core
# - Plot PAE heatmap (rank_1) and save PNG
# - Save a one-row CSV with the metrics to <input_dir_parent>/summary/<name>_pae_summary.csv

import os
import re
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def ask(msg, default=None, allow_empty=False):
    tip = f" [{default}]" if default else ""
    while True:
        s = input(f"{msg}{tip}: ").strip()
        if not s and default:
            return default
        if s or allow_empty:
            return s
        print("Please enter a value.")

def parse_loop_ranges(s):
    # e.g. "392-402;450-460"
    s = s.strip()
    if not s:
        return []
    out = []
    for part in [p.strip() for p in s.split(";") if p.strip()]:
        m = re.match(r"^\s*(-?\d+)\s*[-:]\s*(-?\d+)\s*$", part)
        if not m:
            print(f"Skip unrecognized range: {part}")
            continue
        a, b = int(m.group(1)), int(m.group(2))
        if a > b: a, b = b, a
        out.append((a, b))
    return out

def find_pae_json(colab_dir):
    # collect candidate JSONs
    cands = []
    for base, _, files in os.walk(colab_dir):
        for fn in files:
            low = fn.lower()
            if low.endswith(".json") and "predicted_aligned_error" in low:
                cands.append(os.path.join(base, fn))
            # some colabfold runs store pae as pae.json or similar:
            elif low.endswith(".json") and ("pae" in low) and ("rank" in low or "model" in low):
                cands.append(os.path.join(base, fn))
    if not cands:
        # final fallback: any json that has predicted_aligned_error key
        for base, _, files in os.walk(colab_dir):
            for fn in files:
                if fn.lower().endswith(".json"):
                    cands.append(os.path.join(base, fn))
    if not cands:
        return None

    # prefer names that look like rank_001 / rank_1
    def rank_score(path):
        name = os.path.basename(path).lower()
        score = 0
        if "rank_001" in name or "rank1" in name or "rank_1" in name:
            score += 10
        if "predicted_aligned_error" in name:
            score += 5
        return score

    cands.sort(key=lambda p: rank_score(p), reverse=True)
    # verify JSON has the expected key
    for p in cands:
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                data = json.load(f)
            if "predicted_aligned_error" in data:
                return p, data
        except Exception:
            continue
    return None

def mask_from_ranges(N, ranges):
    mask = np.zeros(N, dtype=bool)
    for a, b in ranges:
        # residue indices are 1-based typically; adjust if needed
        a0 = max(1, a); b0 = min(N, b)
        mask[a0-1:b0] = True
    return mask

def mean_from_mask_pairs(M, maskA, maskB, symmetric=True):
    # M: NxN PAE matrix
    # if symmetric=True, average M[A,B] and M[B,A] together
    if not maskA.any() or not maskB.any():
        return None
    sub1 = M[np.ix_(maskA, maskB)]
    if symmetric:
        sub2 = M[np.ix_(maskB, maskA)]
        vals = np.concatenate([sub1.ravel(), sub2.ravel()])
    else:
        vals = sub1.ravel()
    if vals.size == 0:
        return None
    return float(np.mean(vals))

def plot_pae(M, vmax, title, save_path):
    fig, ax = plt.subplots(figsize=(6,5))
    im = ax.imshow(M, vmin=0, vmax=vmax)  # use matplotlib default colormap
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Predicted alignment error (Å)")
    ax.set_xlabel("Residue index")
    ax.set_ylabel("Residue index")
    ax.set_title(title)
    plt.tight_layout()
    fig.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

# -------- run --------
print("🔧 Minimal PAE Analysis (save to sibling 'summary/')")

colab_dir = ask("ColabFold output directory (absolute path)", "/path/to/colabfold_outputs")
if not os.path.isdir(colab_dir):
    sys.exit(f"Directory not found: {colab_dir}")

name = ask("Name (used for Enzyme & filename prefix)", "BfIMTD")
loop_s = ask("Loop ranges, e.g. 392-402;450-460 (Enter to skip)", "", allow_empty=True)
loop_ranges = parse_loop_ranges(loop_s)

found = find_pae_json(colab_dir)
if not found:
    sys.exit("No suitable PAE JSON found (predicted_aligned_error*.json).")
pae_json_path, data = found

# load matrix
try:
    M = np.array(data["predicted_aligned_error"], dtype=float)
except Exception as e:
    sys.exit(f"Failed to parse predicted_aligned_error matrix: {e}")
vmax = float(data.get("max_predicted_aligned_error", 30.0))
N = M.shape[0]

# define masks
loop_mask = mask_from_ranges(N, loop_ranges) if loop_ranges else np.zeros(N, dtype=bool)
core_mask = ~loop_mask

# compute metrics
mean_global = float(np.mean(M))
mean_ll = mean_from_mask_pairs(M, loop_mask, loop_mask, symmetric=True) if loop_ranges else None
mean_lc = mean_from_mask_pairs(M, loop_mask, core_mask, symmetric=True) if loop_ranges else None
mean_cc = mean_from_mask_pairs(M, core_mask, core_mask, symmetric=True) if loop_ranges else float(np.mean(M))

# save outputs to sibling "summary"
summary_dir = os.path.join(os.path.dirname(colab_dir.rstrip(os.sep)), "summary")
os.makedirs(summary_dir, exist_ok=True)

# PNG
png_path = os.path.join(summary_dir, f"{name}_rank1_PAE_heatmap.png")
plot_pae(M, vmax=vmax, title=f"{name} PAE (rank_1)", save_path=png_path)

# CSV (one row)
row = {
    "Enzyme": name,
    "Residues_N": N,
    "PAE_global_mean_A": round(mean_global, 2),
    "PAE_loop_loop_mean_A": None if mean_ll is None else round(mean_ll, 2),
    "PAE_loop_core_mean_A": None if mean_lc is None else round(mean_lc, 2),
    "PAE_core_core_mean_A": None if mean_cc is None else round(mean_cc, 2),
    "PAE_json_file": os.path.basename(pae_json_path),
    "PAE_png_file": os.path.basename(png_path)
}
df = pd.DataFrame([row])
csv_path = os.path.join(summary_dir, f"{name}_pae_summary.csv")
df.to_csv(csv_path, index=False)

print("✅ Saved:")
print(" -", png_path)
print(" -", csv_path)