#!/usr/bin/env bash

set -e  # if any command fails, quit

REPO="${@:-all}"
USERNAME=$(<.username)

if [[ $REPO == "all" ]]; then
   REPOS=("autochem" "autoio" "autofile" "mechanalyzer")
else
   REPOS=(${REPO})
fi

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

for repo in ${REPOS[@]}; do
   echo cmd: git subrepo status src/_${repo}
   git subrepo status src/_${repo} | while read line; do
      if [[ ${line} == "Remote URL:"*"Auto-Mech/"* || ${line} == "Tracking Branch: dev" ]]; then
         echo -e "${GREEN}${line}${NC}"
      elif [[ ${line} == "Remote URL:"* || ${line} == "Tracking Branch"* ]]; then
         echo -e "${RED}${line}${NC}"
      else
         echo -e "${line}"
      fi
   done
done
