# ==== 用户设置区域 ====
select loop_residues, resi 392-402
select ligand_atoms, ligand_docked
# ======================

delete loop_ligand_contacts

# 选出 loop 原子
select loop_atoms, loop_residues

# 显示 cartoon 和 sticks
show cartoon, all
color slate, all
show sticks, ligand_atoms
color orange, ligand_atoms
show sticks, loop_atoms

dist hbonds, loop_atoms, ligand_atoms, 3.5, 60, 2

# 输出并高亮靠近 ligand 的 loop 残基
iterate (loop_atoms within 6.0 of ligand_atoms), print(resi)
select contact_loop_residues, loop_atoms within 6.0 of ligand_atoms
show sticks, contact_loop_residues
color red, contact_loop_residues

# 美化可视化
set dash_width, 2
set dash_color, yellow
zoom loop_atoms or ligand_atoms