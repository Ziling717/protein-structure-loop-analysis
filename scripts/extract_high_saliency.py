import os
import re
import csv

import os

input_dir = "/Volumes/APFS/protein_prediction/data/deepfri_viz/viz"
print("✅ 路径是否存在：", os.path.exists(input_dir))
print("📄 路径内容：", os.listdir(input_dir) if os.path.exists(input_dir) else "路径不可访问")

# ========== MODIFY THIS ==========
input_dir = "/Volumes/APFS/protein_prediction/data/deepfri_viz/viz"  # 替换为你自己的文件夹路径
output_csv = "saliency_summary.csv"
threshold = 0.8
# =================================

all_data = []

# 遍历所有viz脚本
for filename in os.listdir(input_dir):
    if filename.startswith("pymol_viz_") and filename.endswith(".py"):
        filepath = os.path.join(input_dir, filename)

        # 解析蛋白ID和GO term
        match = re.match(r"pymol_viz_(.+?)_GO_(.+?)\.py", filename)
        if not match:
            continue
        protein_id = match.group(1)
        go_term = "GO:" + match.group(2).replace("_", ":")

        with open(filepath, "r") as f:
            contents = f.read()

        # 提取 stored.cam 列表
        cam_match = re.search(r"stored\.cam\s*=\s*\[([^\]]+)\]", contents)
        if not cam_match:
            continue
        cam_values = cam_match.group(1).strip().split(",")

        # 转换为 float 并提取大于阈值的残基编号
        for i, val in enumerate(cam_values):
            try:
                score = float(val)
                if score > threshold:
                    all_data.append([protein_id, go_term, i + 1, score])  # +1 for 1-based indexing
            except:
                continue

# 写入CSV文件
with open(output_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["protein_id", "go_term", "residue_index", "saliency_score"])
    # 排序：先按 protein_id, 再按 go_term, 再按 saliency_score 降序
    all_data.sort(key=lambda x: (x[0], x[1], -x[3]))
    writer.writerows(all_data)

print(f"✅ Done. Extracted data saved to: {output_csv}")