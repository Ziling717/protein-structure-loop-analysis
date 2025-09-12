# draw_axes_with_labels.pml

# origin
pseudoatom origin, pos=[0,0,0]

# X axis（red）+ label
pseudoatom x_axis, pos=[5,0,0]
distance x_line, origin, x_axis
label x_axis, "X"
color red, x_line
set label_color, red, x_axis

# Y axis（green）+ label
pseudoatom y_axis, pos=[0,5,0]
distance y_line, origin, y_axis
label y_axis, "Y"
color green, y_line
set label_color, green, y_axis

# Z axis（blue）+ label
pseudoatom z_axis, pos=[0,0,5]
distance z_line, origin, z_axis
label z_axis, "Z"
color blue, z_line
set label_color, blue, z_axis

# clear atom balls
hide nonbonded, origin x_axis y_axis z_axis

# make lines thicker
set dash_width, 2.0