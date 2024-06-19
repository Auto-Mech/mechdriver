# MechDriver

## Description

This repository houses the main driver for executing an AutoMech workflow.

## Installation

### Users

End users can install the code using Conda, Mamba, or Pixi.
```
conda install automech # option 1
mamba install automech # option 2
pixi add automech      # option 3
```
Before running the above command, you will need to add `auto-mech` to your list of
channels:
```
conda config --append channels auto-mech # option 1 and 2
pixi project channel add auto-mech       # option 3
```
If `conda-forge` isn't the default channel for your Conda/Mamba installation, you will
also need to set this additional configuration using the command above with `--prepend`.

### Developers

Developers who wish to contribute to the code should follow the instructions
[here](https://github.com/avcopan/amech-dev?tab=readme-ov-file#automech-developer-set-up)
to get set up.
