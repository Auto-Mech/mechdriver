#!/usr/bin/env bash

set -e  # if any command fails, quit

REPO="${@:-all}"
USERNAME=$(<.username)

if [[ $REPO == "all" ]]; then
   REPOS=("autochem" "autoio" "autofile" "mechanalyzer")
else
   REPOS=(${REPO})
fi

for repo in ${REPOS[@]}; do
   # 1. Sync to make sure the fork is up-to-date
   echo cmd: gh repo sync ${USERNAME}/${repo}
   gh repo sync ${USERNAME}/${repo}
   echo

   # 2. Pull from upstream first
   echo cmd: git subrepo pull src/_${repo} -r https://github.com/Auto-Mech/${repo}.git
   git subrepo pull src/_${repo} -r https://github.com/Auto-Mech/${repo}.git
   echo

   # 3. Pull from fork, in case it is ahead of upstream
   echo cmd: git subrepo pull src/_${repo} -r git@github.com:${USERNAME}/${repo}.git
   git subrepo pull src/_${repo} -r git@github.com:${USERNAME}/${repo}.git
   echo
done
