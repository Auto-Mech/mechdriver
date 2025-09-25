#!/usr/bin/env bash

set -e  # if any command fails, quit

if [[ -z "$1" ]]; then
    echo "Must specify 'all' or individual repo to pull from."
    exit 1
fi

REPO=${1}
shift 1
ARGS="$*"  # Additional arguments for git subrepo pull
USERNAME=$(git config --global user.name)

if [[ $REPO == "all" ]]; then
   REPOS=("autochem" "autoio" "autofile" "mechanalyzer")
else
   REPOS=(${REPO})
fi

# Check GitHub CLI authentication
if [[ -z ${ARGS} ]]; then
   if gh auth status; then
      echo "You are already logged into the GitHub CLI. Make sure you have workflow permissions."
      echo "If not, you can run 'pixi run gh-auth' to add them."
   else
      echo "Logging into GitHub CLI... Please follow the prompts to authenticate."
      gh auth login -s workflow
   fi
fi

for repo in ${REPOS[@]}; do
   if [[ -z ${ARGS} ]]; then
      # 1. Sync to make sure the fork is up-to-date
      echo \$ gh repo sync ${USERNAME}/${repo}
      gh repo sync ${USERNAME}/${repo}
      echo

      # 2. Pull from upstream first
      echo \$ git subrepo pull src/_${repo} -r https://github.com/Auto-Mech/${repo}.git
      git subrepo pull src/_${repo} -r https://github.com/Auto-Mech/${repo}.git
      echo
   fi

   # 3. Pull from fork, in case it is ahead of upstream
   echo \$ git subrepo pull src/_${repo} -r git@github.com:${USERNAME}/${repo}.git ${ARGS}
   git subrepo pull src/_${repo} -r git@github.com:${USERNAME}/${repo}.git ${ARGS}
   echo
done
