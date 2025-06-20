## 0. comment block
# ${comment}
## 1. memory block
mem=${memory}GB
## 2. method block
calc=${corr_method}
basis=${basis}
scftype=${scf_method}
## job options block
% if job_options:
${job_options}
% endif
## n. scf options block
% if scf_options:
${scf_options}
% endif
## n. corr options block
%if corr_options:
${corr_options}
% endif
## x. molecule block
charge=${charge}
mult=${mult}
## x. job block
% if job_key == 'optimization':
gopt=full
% elif job_key == 'hessian':
freq=on
% endif

geom=${coord_sys}
${geom}

% if zmat_val_str:
${zmat_val_str}
% endif
