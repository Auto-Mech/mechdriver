## Set the number of atoms 
Number_of_Atoms:                ${natoms}
## Keywords that will generally not change, except for IRCs
Act_energy(kcal/mol):           0.0
Initial_Temperature:            50
Temperature_steps:              40
Temperature_increment:          40
Delta_Energy_rea:               0.
Delta_Energy_pro:               0.
Maxstep:                        ${nsteps}
Npointsint:                     5
Maxtdev:                        0.5
Rearrange(1=yes,0=no)           1
SaddlePoint                     ${saddle_idx}
## Set if projections to be done in internal or cartesian coordinates
% if coord_proj == 'internal':
internalcoord(1=yes)            1
% endif
ds(1=noexp,0=standard)          0
isct_vtst(1=vtst_sct,0=sct)     1
zerocurvature(1)                0
reduced_mass                    1.0
minimum_frequency               50
anim_freq(if_Maxstep=1)         2
## Set projection of the reaction coordinate
% if proj_rxn_coord:
onlyrotors(0=yes,1=no)          1
proj_rea_coo(0=yes(def),1=no)   0
% else:
onlyrotors(0=yes,1=no)          0
proj_rea_coo(0=yes(def),1=no)   1
% endif
## Define all of the rotors
numrotors                       ${nrotors}
% if nrotors > 0:
${rotors_str}
% endif
## Write the geometries, gradients, and Hessians for all Steps
${data_str}\
