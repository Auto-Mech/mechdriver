#!/usr/bin/env bash

set -e  # if any command fails, quit

if [[ -z "$1" ]]; then
    echo "Must specify 'all' or individual repo to pull from."
    exit 1
fi

REPO=${1}
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
