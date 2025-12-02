#!/usr/bin/env bash

# This script updates each repo against a remote
# (It also updates the amech-dev repo.)
#
# Performs a pull --rebase
#
# Arguments:
#   - Remote to update against (default: upstream)
#   - Branch to update (default: dev)

set -e  # if any command fails, quit
REPOS=("autochem" "autoio" "autofile" "mechanalyzer" "mechdriver")

# 0. Read arguments
REMOTE=${1:-upstream}
BRANCH=${2:-dev}

echo "The following commands will be run in each repository:"
echo "    git checkout ${BRANCH}"
echo "    git pull --rebase ${REMOTE} ${BRANCH}"
read -p "Is this what you want to do? [y/n] " yn

if [[ $yn =~ ^[Yy]$ ]]; then
    # 1. Navigate to mechdriver parent directory
    (
        cd ..

        # 2. Loop through each repo and update
        for repo in ${REPOS[@]}
        do
            printf "\n*** Updating in $(realpath ${repo}) ***\n"
            (
                cd ${repo} && \
                git checkout ${BRANCH} && \
                git pull --rebase ${REMOTE} ${BRANCH}
            )
            printf "******\n"
        done
    )
fi