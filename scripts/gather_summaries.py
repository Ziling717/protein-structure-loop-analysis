#!/usr/bin/env python3
import os, glob, pandas as pd

def read_csv_loose(path):
    for enc in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            return pd.read_csv(path, encoding=enc)
        except Exception:
            continue

    return pd.read_csv(path, encoding="utf-8", on_bad_lines="skip")

def main(summary_dir):
    os.makedirs(summary_dir, exist_ok=True)

    # 1) summary pLDDT
    plddt_paths = sorted(glob.glob(os.path.join(summary_dir, "*plddt*summary*.csv")))
    if plddt_paths:
        frames = []
        for p in plddt_paths:
            df = read_csv_loose(p)
            df["__source__"] = os.path.basename(p)
            frames.append(df)
        pd.concat(frames, ignore_index=True)\
          .to_csv(os.path.join(summary_dir, "plddt_all_models.csv"), index=False)
        print("✅ saved:", os.path.join(summary_dir, "plddt_all_models.csv"))
    else:
        print("⚠️ no pLDDT csv found in", summary_dir)

    # 2) summary PAE
    pae_paths = sorted(glob.glob(os.path.join(summary_dir, "*pae*summary*.csv")))
    if pae_paths:
        frames = []
        for p in pae_paths:
            df = read_csv_loose(p)
            df["__source__"] = os.path.basename(p)
            frames.append(df)
        pd.concat(frames, ignore_index=True)\
          .to_csv(os.path.join(summary_dir, "pae_all.csv"), index=False)
        print("✅ saved:", os.path.join(summary_dir, "pae_all.csv"))
    else:
        print("⚠️ no PAE csv found in", summary_dir)

if __name__ == "__main__":

    summary_dir = input("Path to summary/ directory [./summary]: ").strip() or "./summary"
    main(summary_dir)