# use target atoms to create a pocket box
# visualize the box in pymol

# === Step 1: set box center and size ===
center_x = -1.753999948501587
center_y = 6.751999855041504
center_z = -14.369999885559082

size_x = 20
size_y = 30
size_z = 18

# === Step 2: caliculate box boundary ===
x1 = center_x - size_x/2
x2 = center_x + size_x/2
y1 = center_y - size_y/2
y2 = center_y + size_y/2
z1 = center_z - size_z/2
z2 = center_z + size_z/2

# === Step 3: create 8 corner point ===
cmd.pseudoatom("p1", pos=[x1,y1,z1])
cmd.pseudoatom("p2", pos=[x2,y1,z1])
cmd.pseudoatom("p3", pos=[x2,y2,z1])
cmd.pseudoatom("p4", pos=[x1,y2,z1])
cmd.pseudoatom("p5", pos=[x1,y1,z2])
cmd.pseudoatom("p6", pos=[x2,y1,z2])
cmd.pseudoatom("p7", pos=[x2,y2,z2])
cmd.pseudoatom("p8", pos=[x1,y2,z2])

# === Step 4: connect the points to form lines ===
cmd.distance("line1", "p1", "p2")
cmd.distance("line2", "p2", "p3")
cmd.distance("line3", "p3", "p4")
cmd.distance("line4", "p4", "p1")

cmd.distance("line5", "p5", "p6")
cmd.distance("line6", "p6", "p7")
cmd.distance("line7", "p7", "p8")
cmd.distance("line8", "p8", "p5")

cmd.distance("line9",  "p1", "p5")
cmd.distance("line10", "p2", "p6")
cmd.distance("line11", "p3", "p7")
cmd.distance("line12", "p4", "p8")

# === Step 5: clear relunctant lines ===
cmd.hide("nonbonded", "p*")
cmd.set("dash_width", 2.0)
cmd.set("dash_color", "yellow")