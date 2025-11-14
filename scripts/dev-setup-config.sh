#!/usr/bin/env bash

# This script sets up configuration for development.

set -e  # if any command fails, quit

echo
echo Turning off Python output buffering by adding this line to $HOME/.bashrc:
echo "export PYTHONUNBUFFERED=1"
echo "export PYTHONUNBUFFERED=1" >> $HOME/.bashrc
