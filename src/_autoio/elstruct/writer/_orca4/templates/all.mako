## 0. comment block
# ${comment}
## 1. machine options block
%% pal nprocs ${nprocs} end
%% MaxCore ${memory} 
## 2. theoretical method block
% if reference:
! ${reference} ${method} ${basis}
% else:
! ${method} ${basis}
% endif
## 2. job options block (when we get to it)
% if job_key == 'optimization':
! Opt
% elif job_key == 'gradient':
! EnGrad
% elif job_key == 'hessian' and not numerical:
! AnFreq
% elif job_key == 'hessian' and numerical:
! NumFreq
% endif
## 3. theory options block
% if scf_options:
%% scf
${scf_options}
end
% endif
## 4. molecule block
* ${coord_sys} ${charge} ${mult}
${geom}
*

