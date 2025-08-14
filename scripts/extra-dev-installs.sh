#!/usr/bin/env bash

PREFIX=$1

if [ -z "${1}" ]; then
    read -p "Install prefix (default: $HOME/bin): " PREFIX
    PREFIX=${PREFIX:-$HOME/bin}
else
    echo Install prefix: $PREFIX
fi

# # 1. Install HyperQueue
HQ_TGZ_URL="https://github.com/It4innovations/hyperqueue/releases/download/v0.20.0/hq-v0.20.0-linux-x64.tar.gz"
echo
echo Downloading HyperQueue executable into $PREFIX from $HQ_TGZ_URL
curl -sL $HQ_TGZ_URL | tar xzf - -C $PREFIX

echo
echo Putting HyperQueue in your path by adding this line to $HOME/.bashrc:
echo "export PATH=\$PATH:$PREFIX"
echo "export PATH=\$PATH:$PREFIX" >> ~/.bashrc

# 2. Install Git Subrepo
GSR_REPO_URL="https://github.com/ingydotnet/git-subrepo"
echo
echo Downloading Git Subrepo repository from $GSR_REPO_URL
git clone $GSR_REPO_URL $PREFIX/git-subrepo
if [ $? -eq 0 ]; then
    # command succeeded
    echo
    echo Putting Git Subrepo in your path by adding this line to $HOME/.bashrc:
    echo "source $PREFIX/git-subrepo/.rc"
    echo "source $PREFIX/git-subrepo/.rc" >> ~/.bashrc
else
    # command failed
    echo
    echo Git Subrepo appears to be installed. Make sure you have this line in $HOME.bashrc:
    echo "source $PREFIX/git-subrepo/.rc"
fi
