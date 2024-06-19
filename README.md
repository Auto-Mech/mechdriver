# MechDriver

## Description

This repository houses the main driver for executing an AutoMech workflow.

## Installation

### Developers

The AutoMech code is in development, so we encourage users to install it in developer mode and send us bug reports.
See
[here](https://github.com/avcopan/amech-dev?tab=readme-ov-file#automech-developer-set-up)
for instructions to get set up.


### Users

End-users who are unable to contribute (see above) can install the code using Conda, Mamba, or Pixi.
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
