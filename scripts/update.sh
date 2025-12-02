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
REPOS=("autochem" "autoio" "autofile" "mechanalyzer")

# 1. Loop through each repo and update
for repo in ${REPOS[@]}
do
    version=$(pixi run --manifest-path ../$repo current-version)
    echo Setting $repo version to $version
    sed -i -E "s/($repo *= *\"==)[0-9]+\.[0-9]+\.[0-9]+(\" *)/\1${version}\2/" pyproject.toml
done

# 2. Update lockfile
echo Updating lockfile
pixi lock