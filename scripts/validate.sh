#!/usr/bin/env bash

set -e  # if any command fails, quit

REPO="${@:-all}"

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

EXIT_CODE=0
for repo in ${REPOS[@]}; do
   echo Validating ${repo}...
   remote=$(cat src/_${repo}/.gitrepo | sed -n -E -e 's/^\s*remote = (.*)/\1/p')
   branch=$(cat src/_${repo}/.gitrepo | sed -n -E -e 's/^\s*branch = (.*)/\1/p')
   commit=$(cat src/_${repo}/.gitrepo | sed -n -E -e 's/^\s*commit = (.*)/\1/p')
   ping=$(git ls-remote https://github.com/Auto-Mech/${repo}.git | sed -n -E -e "s/^(${commit})/\1/p")
   if [[ ${remote} == *"Auto-Mech/"* && ${branch} == "dev" ]]; then
      if [[ -n "${ping}" ]]; then
         echo -e "${GREEN}  OK: ${remote} ${branch} ${commit}${NC}"
      else
         echo -e "${RED}  PING FAILED: ${remote} ${branch} ${commit}${NC}"
         EXIT_CODE=1
      fi
   else
      echo -e "${RED}  NOT SYNCED: ${remote} ${branch} ${commit}${NC}"
      EXIT_CODE=1
   fi
done

if [[ ${EXIT_CODE} -eq 0 ]]; then
   echo
   echo PASSED
else
   echo
   echo FAILED
fi

exit ${EXIT_CODE}
