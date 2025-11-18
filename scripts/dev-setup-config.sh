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

# 2. Add mechenv function to .bashrc
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
MECHDRIVER_PATH="$( realpath "$SCRIPT_DIR/.." )"

ENV_FUNCTION_BLOCK=$(cat <<EOF
mechenv() {
    local env_name="\${1:-dev}"
    eval "\$(pixi shell-hook -e "\$env_name" --manifest-path "$MECHDRIVER_PATH")"
}
EOF
)

echo
if ! grep -Fq "mechenv()" "$HOME/.bashrc"; then
    echo "Adding mechenv() function to $HOME/.bashrc:"
    echo "$ENV_FUNCTION_BLOCK"
    echo "$ENV_FUNCTION_BLOCK" >> "$HOME/.bashrc"
else
    echo "Function mechenv() is already present in $HOME/.bashrc"
fi
