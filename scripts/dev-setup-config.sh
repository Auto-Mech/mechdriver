#!/usr/bin/env bash

# This script sets up configuration for development.

set -e  # if any command fails, quit

# 1. Turn off Python output buffering
BUFFER_CONFIG_LINE='export PYTHONUNBUFFERED=1'
echo Turning off Python output buffering by adding this line to $HOME/.bashrc:
echo $BUFFER_CONFIG_LINE
echo $BUFFER_CONFIG_LINE >> $HOME/.bashrc

# 2. Add alias to activate mechdriver environment
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MECHDRIVER_PATH="$( realpath "$SCRIPT_DIR/.." )"
ENV_CONFIG_LINE="alias mechenv='eval \"\$(pixi shell-hook -e dev --manifest-path $MECHDRIVER_PATH)\"'"
echo
echo Adding alias 'mechenv' for activating dev environment to $HOME/.bashrc:
echo $ENV_CONFIG_LINE
echo $ENV_CONFIG_LINE >> $HOME/.bashrc
