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
    echo cmd: git subrepo push src/_${repo} -r git@github.com:${USERNAME}/${repo}.git
    git subrepo push src/_${repo} -r git@github.com:${USERNAME}/${repo}.git
done
