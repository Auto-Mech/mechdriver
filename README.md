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

## Running

### Basic

To run the main AutoMech workflow, you will need an `inp/` directory with input files as in one of the examples [here](https://github.com/Auto-Mech/amech-dev/tree/main/examples). You can then run the main workflow as follows:
```
automech run &> out.log &
```
To see live output from the workflow, you will need to turn off Python buffering before you run the above command.
```
export PYTHONUNBUFFERED=1
```
You can then refresh your shell before running the above command and run `tail -f out.log` to see the live output and monitor progress.

### Subtasks

Workflow parallelization is currently not automated in AutoMech. However, if you are on a cluster with direct SSH node access and permissions to run, you can run the following commands to split an AutoMech workflow into subtasks and run them in parallel.

(1.) You can set-up these subtask jobs as follows:
```
automech subtasks setup
```
This will parse your `inp/` directory and create individual subdirectories for running each individual task for each individual species or reaction/TS. These directories will go in a folder called `subtasks/`.

(2.) If you are using the [amech-dev](https://github.com/Auto-Mech/amech-dev) Pixi environment, you can run the subtasks in parallel on a list of nodes as follows:
```
pixi run subtasks csed-0008 csed-0009 csed-0010  # or csed-00{08..10}
```
Otherwise, you can execute the subtasks in parallel as follows:
```
automech subtasks run-adhoc -n <comma-separated list of nodes> -a <environment activation commands>
```
Here, the list of nodes might be `csed-0008,csed-0009,csed-0010` and the activation command might be `"$(pixi shell-hook)"`.

(3.) To check the progress of your subtask run, you can use the following command:
```
automech subtasks status
```
This will print a color-coded table showing which tasks have failed for which species/reactions. It will also generate a `check.log` file with the paths to log files that have have not completed successfully or have a warning.

