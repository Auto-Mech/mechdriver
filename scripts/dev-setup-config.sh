#!/usr/bin/env bash

# This script sets up configuration for development.

set -e  # if any command fails, quit

# 1. Turn off Python output buffering
BUFFER_CONFIG_LINE='export PYTHONUNBUFFERED=1'
if ! grep -Fxq "$BUFFER_CONFIG_LINE" $HOME/.bashrc; then
    echo Turning off Python output buffering by adding this line to $HOME/.bashrc:
    echo $BUFFER_CONFIG_LINE
    echo $BUFFER_CONFIG_LINE >> $HOME/.bashrc
else
    echo Python output buffering is already turned off in $HOME/.bashrc
fi

# 2. Add alias to activate mechdriver environment
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MECHDRIVER_PATH="$( realpath "$SCRIPT_DIR/.." )"
ENV_CONFIG_LINE="alias mechenv='eval \"\$(pixi shell-hook -e dev --manifest-path $MECHDRIVER_PATH)\"'"
echo
if ! grep -Fxq "$ENV_CONFIG_LINE" $HOME/.bashrc; then
    echo Adding alias 'mechenv' for activating dev environment to $HOME/.bashrc:
    echo $ENV_CONFIG_LINE
    echo $ENV_CONFIG_LINE >> $HOME/.bashrc
else
    echo Alias 'mechenv' is already present in $HOME/.bashrc
fi
