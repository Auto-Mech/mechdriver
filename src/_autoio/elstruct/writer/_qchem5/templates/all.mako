$comment
${comment}
$end

$molecule
${charge} ${mult}
${geom}

% if zmat_var_vals != '':
${zmat_var_vals}
% endif
$end

% if zmat_const_vals != '':
$opt
CONSTRAINT
${zmat_const_vals}
ENDCONSTRAINT
$end
% endif

$rem
    JOBTYPE           ${job_key}
% if job_key == 'optimization':
    GEOM_OPT_COORDS   2
% endif
    METHOD            ${method} 
    UNRESTRICTED      ${unrestricted}
    BASIS             ${basis}
% if xcgrid is not None:
    XC_GRID           ${xcgrid}
% endif
    MEM_TOTAL         ${memory}
$end
