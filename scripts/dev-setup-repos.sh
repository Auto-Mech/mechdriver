#!/usr/bin/env bash

# This script downloads each forked repo from the user's GitHub, checks out the dev
# branch and updates it against the main Auto-Mech repo
#
# Arguments:
#   - GitHub username (required)
#   - Repo branch (default: dev)

set -e  # if any command fails, quit
REPOS=("autochem" "autoio" "autofile" "mechanalyzer" "mechdriver")

echo "Setting up GitHub repositories: ${REPOS[@]}"

# 0. Read arguments
USERNAME=${1}
UPSTREAM=${2}
MODE=${3}
BRANCH=${4:-dev}

DEFAULT_USERNAME=$(git config --global user.name)
if [[ -z "$DEFAULT_USERNAME" && -z "$USERNAME" ]]; then
    read -p "Please enter your GitHub username: " USERNAME
fi

echo "Press enter to choose the default values..."
if [ -z "$USERNAME" ]; then
    read -p "  Git username (${DEFAULT_USERNAME} [default] or enter alternative): " INPUT
    USERNAME=${INPUT:-$DEFAULT_USERNAME}
fi
if [ -z "${2}" ]; then
    read -p "  Update against Auto-Mech upstream? (yes [default] or no): " UPSTREAM
    UPSTREAM=${UPSTREAM:-yes}
fi
if [ -z "${3}" ]; then
    read -p "  How would you like to clone? (ssh [default] or http): " MODE
    MODE=${MODE:-ssh}
fi

echo "Arguments:"
echo "  Username - ${USERNAME}"
echo "  Update   - ${UPSTREAM}"
echo "  Mode     - ${MODE}"
echo "  Branch   - ${BRANCH}"
if [[ -n "$1" && -n "$2" && -n "$3" ]]; then
    echo "All arguments were provided via command line. Skipping confirmation."
else
    read -p "Is this correct? If so, press enter to continue"
fi

CLONE_PREFIX="https://github.com/${USERNAME}"
if [ "${MODE}" == "ssh" ]; then
    CLONE_PREFIX="git@github.com:${USERNAME}"
fi

# 1. Navigate to mechdriver parent directory
(
    cd ..

    # 2. Loop through each repo and download it
    for repo in ${REPOS[@]}
    do
        printf "\n*** Cloning from ${CLONE_PREFIX}/${repo}.git\n"
        if [ -d "${repo}" ]; then
            # a. If the directory already exists, skip it
            echo Directory ${repo} already exists. Skipping...
        else
            # a. Clone the repo
            git clone ${CLONE_PREFIX}/${repo}.git
            # b. If it worked, enter the repo, add Auto-Mech as a remote, and add the branch
            # both locally and on GitHub
        fi
        if [ "${UPSTREAM}" == "yes" ]; then
            (
                # i. Enter the repository
                cd ${repo}
                # ii. If the desired branch isn't the default one, fetch it from origin and
                # switch to it
                branch=$( git branch | tr -d [*] | xargs )
                if [[ ${branch} != ${BRANCH} ]]; then
                    git fetch origin ${BRANCH} && \
                    git branch ${BRANCH} FETCH_HEAD && \
                    git checkout ${BRANCH}
                fi
                # iii. Rebase the selected branch against upstream
                git remote add upstream https://github.com/Auto-Mech/${repo} || true
                git pull --rebase upstream ${BRANCH}
            )
        fi
        printf "***\n"
    done
)
