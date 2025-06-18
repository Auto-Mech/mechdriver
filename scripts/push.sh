#!/usr/bin/env bash

set -e  # if any command fails, quit

REPO=${1:-"all"}
shift 1
ARGS="$*"  # Additional arguments for git subrepo push
USERNAME=$(<.username)

if [[ $REPO == "all" ]]; then
   REPOS=("autochem" "autoio" "autofile" "mechanalyzer")
else
   REPOS=(${REPO})
fi

for repo in ${REPOS[@]}; do
    echo \$ git subrepo push src/_${repo} -r git@github.com:${USERNAME}/${repo}.git ${ARGS}
    git subrepo push src/_${repo} -r git@github.com:${USERNAME}/${repo}.git ${ARGS}
done
