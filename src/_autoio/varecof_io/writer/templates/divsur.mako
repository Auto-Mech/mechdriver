# Orientation Frame
# origin: 1
# orientation:
#  3
#  ^
#  |y
#  |  x
#  4--->2
# Frame i j k l
#       ^ ^ ^ ^
#       1 2 3 4

# Pivot Points for Fragment 1
PivotPoints     ${npivot1}
Frame ${frame1}

${pivot_xyz_string1}

# Pivot Points for Fragment 2
PivotPoints     ${npivot2}
Frame ${frame2}

${pivot_xyz_string2}

# interpointal distances
Distances

${dist_coords_string}

Conditions ${nconditions}
${conditions_string}

Cycles ${ncycles}

${r1_string}
% if r2_string != '':
${r2_string}
% endif
% if d1_string != '':
${d1_string}
% endif
% if d2_string != '':
${d2_string}
% endif
% if p1_string != '':
${p1_string}
% endif
% if t1_string != '':
${t1_string}
% endif
% if p2_string != '':
${p2_string}
% endif
% if t2_string != '':
${t2_string}
% endif

EndSurface
