"""
temp file for a routine for the pivot points
"""

geom1 = ()
geom2 = ()
rts = 0
aabs1 = 0
aabs3 = 0
babs1 = 0
babs3 = 0
lr_rdists = ()
sr_rdists = ()
d1dists = ()
d2dists = () 
p1angs = ()  
p2angs = ()  
t1angs = ()  
t2angs = ()  


# Get the position of the pivot atoms for fragment 1
# Check the size of the molecule (def good for natom1 = 2)
xyzP1 = calc_pivot_point_xyz(geom1[0], geom1[1], geom1[2],
                             rts, aabs1, babs1)
rho1, angle1, theta1 = calc_pivot_point_internal(xyzP1, geom1[0])

# Get the position of the pivot atoms for fragment 2
# Consider multiple cases of the number of atoms
if not automol.geom.is_atom(geom2):
    if automol.geom.natoms(geom2) >= 3:
        xyzP2 = calc_pivot_point_xyz(geom2[0], geom2[1], geom2[2],
                                     rts, aabs2, -babs3)
        rho2, angle2, theta2 = calc_pivot_point_internal(xyzP2, geom2[0])

# Calculate the number of pivot points
if automol.geom.is_atom(geom1):
    npivot1 = 1
else:
    npivot1 = 2
if automol.geom.is_atom(geom2):
    npivot2 = 1
else:
    npivot2 = 2

# Write the long-range divsur file
lr_divsur_inp_str = varecof_io.writer.divsur_input(
    lr_rdists, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

# Wrute the short-range divsur file
sr_divsur_inp_str = varecof_io.writer.divsur_input(
    sr_rdists, npivot1, npivot2, xyzP1, xyzP2,
    d1dists=d1dists,
    d2dists=d2dists,
    p1angs=p1angs,
    p2angs=p2angs,
    t1angs=t1angs,
    t2angs=t2angs)

