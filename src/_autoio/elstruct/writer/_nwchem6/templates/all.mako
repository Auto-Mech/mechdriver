## 0. comment block
start ${start_title}
title "${comment}"
## 1. memory block
scratch_dir ${scratch_dir}
memory stack ${memory_stack} mb heap ${memory_heap} mb ${memory_global} 10000 mb noverify
## 2. molecule block
${geom}
Q: how to set charge and spin
## 3. method block
basis
  all library ${basis}
end
## 4. job block
% if job_key == 'energy':
task ${method}
% elif job_key == 'gradient':
task ${method} gradient
% elif job_key == 'hessian':
task ${method} frequencies
% elif job_key == 'optimization':
task ${method} optimize
% endif
