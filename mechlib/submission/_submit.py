""" Library of BASH submission scripts for various jobs
"""


PROJROT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "RPHt.exe >& /dev/null"
)
MESSPF = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=10\n"
    "ulimit -c 0\n"
    "/lcrc/project/CMRP/amech/MESS/build/messpf pf.inp pf.out >> stdout.log &> stderr.log"
)
MESSRATE = (
    "#!/usr/bin/env bash\n"
    "export OMP_NUM_THREADS=8\n"
    "ulimit -c 0\n"
    "/lcrc/project/CMRP/amech/MESS/build/mess mess.inp rate.out >> stdout.log &> stderr.log"
)
VARECOF = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/multi tst.inp >& varecof.out"
)
MCFLUX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/mc_flux mc_flux.inp"
)
TSTCHECK = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "/home/ygeorgi/build/rotd/tst_check"
)
THERMP = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "THERMP FILE"
)
PAC99 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "PACC << EOF >& pac99.out\n"
    "FORMULA\n"
    "EOF"
)
DSARRFIT = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "dsarrfit.x_cfg"
)
G09 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "g09 run.inp run.out >> stdout.log &> stderr.log"
)
PSI4 = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "psi4 -i run.inp -o run.out >> stdout.log &> stderr.log"
)
MOLPRO = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro -n 4 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro --mppx -n 4 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MREF = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro -n 8 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)
MOLPRO_MREF_MPPX = (
    "#!/usr/bin/env bash\n"
    "ulimit -c 0\n"
    "molpro --mppx -n 12 run.inp -o run.out "
    "--nouse-logfile --no-xml-output >> "
    "stdout.log &> stderr.log"
)

DEFAULT_SCRIPT_DCT = {
    'projrot': PROJROT,
    'messpf': MESSPF,
    'messrate': MESSRATE,
    'varecof': VARECOF,
    'mcflux': MCFLUX,
    'tstchk': TSTCHECK,
    'thermp': THERMP,
    'pac99': PAC99,
    'dsarrfit': DSARRFIT,
    'gaussian09': G09,
    'molpro2015': MOLPRO,
    'molpro2015_mppx': MOLPRO_MPPX,
    'molpro2015_mr': MOLPRO_MREF,
    'molpro2015_mr_mppx': MOLPRO_MREF_MPPX
}
