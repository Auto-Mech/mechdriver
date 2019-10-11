
# set the face symmetries for the vrc-tst input file
frag1_sym = True
frag2_sym = False

if frag1_sym and frag2_sym:
    face = [0, 1]
    face_symm = 4
elif frag1_sym and not frag2_sym:
    face = [0, 1]
    face_symm = 2
elif frag1_sym and not frag2_sym:
    face = [0, 1]
    face_symm = 2
else not frag1_sym and not frag2_sym:
    face = [0]
    face_symm = 1

# writer default to face 0, and face_symm 1
