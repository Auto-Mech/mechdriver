## 0. machine options block
! ${comment}
memory,${memory},m
% if machine_options:
${machine_options}
% endif
## 1. general lines block 1
% if gen_lines_1 != '':
${gen_lines_1}
% endif
## 2. molecule block
% if mol_options:
${mol_options}
% endif
angstrom
geometry = {
${geom}
}
% if zmat_vals != '':
${zmat_vals}
% endif
## 3. spin and charge block
set,spin=${spin}
set,charge=${charge}
## 4. basis block
basis=${basis}
## 5. general lines block 2 (for fancy wavefunction guessing)
% if gen_lines_2 != '':
${gen_lines_2}
% endif
## 6. method block 
% if scf_method:
{${scf_method},${scf_options}}
% endif
% if ismultiref:
{casscf,${casscf_options}}
% endif
% if corr_method:
{${corr_method},${corr_options}}
% endif

## 7.job block
% if job_key == 'gradient':
{force,${job_options}}
% elif job_key == 'optimization':
{optg,${job_options}}
status
% elif job_key == 'hessian':
{freq,${job_options}}
put,molden,freq.molden
status
% elif job_key == 'energy':
% if 'mp2' in corr_method:
status,${scf_method}-SCF,crash
% else:
status,all,crash
% endif
% endif

${gen_lines_3}
