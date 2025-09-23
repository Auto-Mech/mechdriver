#!/usr/bin/env bash

IS_TRAJ=0

while getopts "t" opt; do
  case $opt in
    t)
      echo "Option -t was provided. Preparing trajectory visualization..."
      IS_TRAJ=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

shift $((OPTIND-1))

WD=${INIT_CWD:-$(pwd)}
ARG=${1}

SCRIPT_HEAD=$(cat << EOF
# %%
from pathlib import Path
import automol

geo_file = Path("${ARG}")
EOF
)
MAIN_BODY=$(cat << EOF
geo = automol.geom.from_xyz_string(geo_file.read_text())
automol.geom.display(geo)
EOF
)
TRAJ_BODY=$(cat << EOF
geos, comments = zip(
    *automol.geom.from_xyz_trajectory_string(geo_file.read_text()),
    strict=True,
)
automol.geom.display_trajectory(geos)
EOF
)

if [ $IS_TRAJ -eq 0 ]; then
    STEM="vis"
    SCRIPT="${SCRIPT_HEAD}\n${MAIN_BODY}"
else
    STEM="vis_traj"
    SCRIPT="${SCRIPT_HEAD}\n${TRAJ_BODY}"
fi

cd ${WD}
printf "${SCRIPT}" > ${STEM}.py

jupytext --to ipynb ${STEM}.py
rm ${STEM}.py
