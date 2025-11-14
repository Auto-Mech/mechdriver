#!/usr/bin/env bash

# This script downloads HyperQueue and installs it in the user's path.
#
# Arguments:
#   - Install prefix (default: $HOME/bin)
#   - HyperQueue URL (default: <URL of .tar.gz from v0.20.0>)

set -e  # if any command fails, quit

echo "Downloading and installing HyperQueue..."

PREFIX=$1
VERSION=$2

if [ -z "${1}" ]; then
    read -p "Install prefix (default: $HOME/bin): " PREFIX
    PREFIX=${PREFIX:-$HOME/bin}
else
    echo Install prefix: $PREFIX
fi

if [ -z "${2}" ]; then
    read -p "HyperQueue version to install (default: 0.20.0): " VERSION
    VERSION=${VERSION:-"0.20.0"}
else
    echo HyperQueue version: $VERSION
fi

# # 1. Install HyperQueue
HQ_TGZ_URL="https://github.com/It4innovations/hyperqueue/releases/download/v${VERSION}/hq-v${VERSION}-linux-x64.tar.gz"
echo
echo Downloading HyperQueue executable into $PREFIX from $HQ_TGZ_URL
curl -sL $HQ_TGZ_URL | tar xzf - -C $PREFIX

echo
EXPORT_PATH_LINE="export PATH=\$PATH:$PREFIX"
if ! grep -Fxq "$EXPORT_PATH_LINE" $HOME/.bashrc; then
    echo Putting HyperQueue in your path by adding this line to $HOME/.bashrc:
    echo "export PATH=\$PATH:$PREFIX"
    echo "export PATH=\$PATH:$PREFIX" >> $HOME/.bashrc
else
    echo HyperQueue path is already set in $HOME/.bashrc
fi
