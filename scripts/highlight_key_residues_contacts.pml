# ===== User Settings =====
set ligand_resname, IG3
select ligand_atoms, resn IG3

# Define residue groups by residue numbers
select catalytic_residues, resi 421+440
select negative_residues, resi 424+426+444+446+470+508
select positive_residues, resi 348+419+462+492+494

# ===== Display Settings =====
hide everything
show cartoon, all
color slate, all

# Ligand
show sticks, ligand_atoms
color orange, ligand_atoms

# Catalytic
show sticks, catalytic_residues
color red, catalytic_residues
dist cat_contacts, catalytic_residues, ligand_atoms, 4.0, 90, 2
set dash_color, red, cat_contacts

# Negative
show sticks, negative_residues
color blue, negative_residues
dist neg_contacts, negative_residues, ligand_atoms, 4.0, 90, 2
set dash_color, blue, neg_contacts

# Positive
show sticks, positive_residues
color purple, positive_residues
dist pos_contacts, positive_residues, ligand_atoms, 4.0, 90, 2
set dash_color, purple, pos_contacts

# Clean view
hide labels
set dash_width, 2
zoom ligand_atoms or catalytic_residues or negative_residues or positive_residues