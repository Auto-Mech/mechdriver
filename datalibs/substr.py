""" Library of strings for bash submission
"""

# ProjRot
PROJROT = ("#!/usr/bin/env bash\n"
           "RPHt.exe")

# MESS
PF = ("#!/usr/bin/env bash\n"
      "messpf pf.inp build.out >> stdout.log &> stderr.log")
RATE = ("#!/usr/bin/env bash\n"
        "mess mess.inp build.out >> stdout.log &> stderr.log")

# VaReCoF
VARECOF = ("#!/usr/bin/env bash\n"
           "/home/ygeorgi/build/rotd/multi ")
MCFLUX = ("#!/usr/bin/env bash\n"
          "/home/ygeorgi/build/rotd/mc_flux ")
CONV_MULTI = ("#!/usr/bin/env bash\n"
              "/home/ygeorgi/build/rotd/mc_flux ")
TST_CHECK = ("#!/usr/bin/env bash\n"
             "/home/ygeorgi/build/rotd/tst_check ")
