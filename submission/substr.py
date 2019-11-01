""" Library of strings for bash submission to external codes
"""

# ProjRot
PROJROT = ("#!/usr/bin/env bash\n"
           "/lcrc/project/CMRP/pacc/ProjRot/build/RPHt.exe >& /dev/null")
#           "RPHt.exe >& /dev/null")

# MESS
MESSPF = ("#!/usr/bin/env bash\n"
          "export OMP_NUM_THREADS=10\n"
          "messpf pf.inp pf.out >> stdout.log &> stderr.log")
MESSRATE = ("#!/usr/bin/env bash\n"
            "export OMP_NUM_THREADS=10\n"
            "mess mess.inp rate.out >> stdout.log &> stderr.log")

# VaReCoF
VARECOF = ("#!/usr/bin/env bash\n"
           "/home/ygeorgi/build/rotd/multi")
MCFLUX = ("#!/usr/bin/env bash\n"
          "/home/ygeorgi/build/rotd/mc_flux")
CONV_MULTI = ("#!/usr/bin/env bash\n"
              "/home/ygeorgi/build/rotd/mc_flux")
TST_CHECK = ("#!/usr/bin/env bash\n"
             "/home/ygeorgi/build/rotd/tst_check")

# Thermo
THERMP = ("#!/usr/bin/env bash\n"
          "THERMP FILE")
PAC99 = ("#!/usr/bin/env bash\n"
         "PACC << EOF\n"
         "FORMULA\n"
         "EOF")
# NASAPOLY = ("#!/usr/bin/env bash\n"
#             "cp ../PF/build.out pf.dat\n"
#             "cp /tcghome/sjklipp/PACC/nasa/new.groups .\n"
#             "python /tcghome/sjklipp/PACC/nasa/makepoly.py"
#             " >> stdout.log &> stderr.log")

# Electronic structure
MOLPRO_PATH_STR = ('/home/sjklipp/bin/molpro')

# dsarrfit (likely not needed, as calling external not done)
DSARRFIT = ("#!/usr/bin/env bash\n"
            "dsarrfit.x_cfg")
