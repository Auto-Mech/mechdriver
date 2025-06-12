tmpr_grid       50      10       1.05     61     Kelvin  # temperature grid
${ener_grid}
${amom_grid}

pot_smp_max     ${nsamp_max} 
pot_smp_min     ${nsamp_min} 
tot_smp_max     100000
tot_smp_min     1000
pot_smp_len     1

flux_rel_err    ${flux_err}       nu
flux_out_mode   ej-resolved
min_atom_dist   1.5     Bohr
priority_exp    0.      nu

pot_name        molpro          # potential to use
pot_type        molpro          #  potential type
pot_file        molpro.inp
pes_size        ${pes_size}

opt_method      none

sampling        multifacet
face            ${faces}
face_symm       ${faces_symm}
ds_inp_file     divsur.inp
ds_out_file     divsur.out

mol_spec_file   structure.inp
save_file       flux.save
flux_file       flux.dat
therm_flux_out  1
therm_flux_name flux.out

raw_smp_file    sampling.raw
is_smp_out      0

work_dir        scratch/node
log_file        molpro.log
run_time        3000000
save_time       3000
flux_debug      0
comm_debug      0
nice            0
