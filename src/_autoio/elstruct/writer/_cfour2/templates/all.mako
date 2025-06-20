## 0. comment block
${comment}
## 1. molecule block
${geom}

% if zmat_var_vals:
${zmat_var_vals}
% endif
## job options block
% if job_key == 'energy':
*CFOUR(GEO_METHOD=SINGLE_POINT
VIBRATION=0
% elif job_key == 'gradient':
*CFOUR(GEO_METHOD=SINGLE_POINT
DERIV_LEVEL=1
VIBRATION=0
% elif job_key == 'hessian' and not numerical:
*CFOUR(VIBRATION=ANALYTIC
% elif job_key == 'hessian' and numerical:
*CFOUR(VIBRATION=FINDIF
% elif job_key == 'optimization' and not saddle:
*CFOUR(GEO_METHOD=NR
VIBRATION=0
% elif job_key == 'optimization' and saddle:
*CFOUR(GEO_METHOD=TS
VIBRATION=0
% endif
## n. theoretical method block
CALC_LEVEL=${method}
REFERENCE=${reference}
BASIS=${basis}
## n. scf options block
% if scf_options:
${scf_options}
% endif
## n. corr options block
%if corr_options:
${corr_options}
% endif
## n. molecule block
CHARGE=${charge}
MULTIPLICITY=${mult}
## n. coord sys block
COORDS=${coord_sys}
UNITS=ANGSTROM
## n. machine options block
MEMORY_SIZE=${memory}
MEM_UNIT=GB)
## n gen lines

% if gen_lines:
${gen_lines}
% endif
