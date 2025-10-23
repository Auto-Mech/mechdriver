#!/usr/bin/env bash

# This script downloads HyperQueue and installs it in the user's path.
#
# Arguments:
#   - Install prefix (default: $HOME/bin)
#   - HyperQueue URL (default: <URL of .tar.gz from v0.20.0>)

set -e  # if any command fails, quit

echo "Downloading and installing HyperQueue..."

PREFIX=$1
HQ_TGZ_URL=$2

if [ -z "${1}" ]; then
    read -p "Install prefix (default: $HOME/bin): " PREFIX
    PREFIX=${PREFIX:-$HOME/bin}
else
    echo Install prefix: $PREFIX
fi

# # 1. Install HyperQueue
HQ_TGZ_URL=${HQ_TGZ_URL:-"https://github.com/It4innovations/hyperqueue/releases/download/v0.20.0/hq-v0.20.0-linux-x64.tar.gz"}
echo
echo Downloading HyperQueue executable into $PREFIX from $HQ_TGZ_URL
curl -sL $HQ_TGZ_URL | tar xzf - -C $PREFIX

echo
echo Putting HyperQueue in your path by adding this line to $HOME/.bashrc:
echo "export PATH=\$PATH:$PREFIX"
echo "export PATH=\$PATH:$PREFIX" >> ~/.bashrc
