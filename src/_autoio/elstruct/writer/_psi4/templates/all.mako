## 0. machine options block
#! ${comment}
memory ${memory} GB
${machine_options}
## 1. molecule block
molecule {
${mol_options}
${charge} ${mult}
${geom}
% if zmat_vals:
${zmat_vals}
% endif
}
% if frozen_dis_strs or frozen_ang_strs or frozen_dih_strs:
set optking {
  frozen_distance = ("
${'\n'.join('    {}'.format(dis_str) for dis_str in frozen_dis_strs)}
  ")
  frozen_bend = ("
${'\n'.join('    {}'.format(ang_str) for ang_str in frozen_ang_strs)}
  ")
  frozen_dihedral = ("
${'\n'.join('    {}'.format(dih_str) for dih_str in frozen_dih_strs)}
  ")
}
% endif
## 2. theory block
set basis ${basis}
set reference ${reference}
set scf_type ${scf_typ}
set mp2_type ${mp2_typ}
% if scf_options:
${scf_options}
% endif
% if corr_options:
${corr_options}
% endif
% if job_options:
${job_options}
% endif
% if gen_lines != '':
${gen_lines}
% endif
## 3. job options block
% if job_key == 'energy':
energy('${method}')
% elif job_key == 'optimization' and method == 'hf' and reference == 'rohf':
optimize('${method}', dertype=1)
% elif job_key == 'optimization' and method == 'hf' and reference == 'uhf':
set {
  intrafrag_step_limit     0.1
  intrafrag_step_limit_min 0.1
  intrafrag_step_limit_max 0.1
}
set full_hess_every 0
optimize('${method}', dertype=1)
% elif job_key == 'optimization':
optimize('${method}')
% elif job_key == 'gradient':
_, wfn = gradient('${method}', return_wfn=True)
grad = wfn.gradient()
grad.name = 'Gradient'
grad.print_out()
% elif job_key == 'hessian' and method == 'hf' and reference == 'rohf':
_, wfn = frequency('${method}', return_wfn=True, dertype=1)
% elif job_key == 'hessian' and method == 'hf' and reference == 'uhf':
_, wfn = frequency('${method}', return_wfn=True, dertype=1)
% elif job_key == 'hessian':
_, wfn = frequency('${method}', return_wfn=True)
% elif job_key == 'irc' and method == 'hf' and reference == 'rohf':
set 
set opt_type irc
set irc_direction ${irc_direction}
_, wfn = frequency('${method}', return_wfn=True, dertype=1)
optimize('${method}', dertype=1)
% elif job_key == 'irc' and method == 'hf' and reference == 'uhf':
set opt_type irc
set irc_direction ${irc_direction}
_, wfn = frequency('${method}', return_wfn=True, dertype=1)
optimize('${method}', dertype=1)
% elif job_key == 'irc':
set opt_type irc
set irc_direction ${irc_direction}
_, wfn = frequency('${method}', return_wfn=True)
optimize('${method}')
% elif job_key == 'molecular_properties':
properties('${method}', properties=['dipole', 'polarizability'])
% endif
% if job_key == 'hessian':
hess = wfn.hessian()
hess.name = 'Hessian'
hess.print_out()
% endif
