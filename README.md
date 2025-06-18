# MechDriver

This repository houses the main driver for executing an AutoMech workflow.
The lower-level submodules of AutoMech are also included in this repository using
[git-subrepo](https://github.com/ingydotnet/git-subrepo), so it can be used to install
the full suite.

## Install

### Developer Mode

If you haven't already, follow the instructions
[here](https://pixi.sh/latest/installation/) to install the Pixi package manager.
It should only take about 30 seconds.
Then turn off Python output buffering by adding the following environment variable to your `~/.bashrc`.
```
export PYTHONUNBUFFERED=1
```

Installing the code in developer mode is then a three step process.

1. Fork this repository and clone the fork to your local machine.
```
git clone git@github.com:<username>/mechdriver.git
cd mechdriver
```
2. In your cloned repository, use Pixi to create the MechDriver developer environment,
which includes all of the necessary dependencies.
```
pixi install -e dev
```
3. Check that the installation worked by running the following help command.
```
pixi shell  # activate the environment
automech --help
```

Once installed, you can activate this environment from anywhere on your system as
follows.
```
pixi shell --manifest-path /path/to/mechdriver
```
If activating within a shell script, the above command wil not work.
In that case, you can activate as follows.
```
eval "$(pixi shell-hook --manifest-path /path/to/mechdriver)"
```

## Run

### Examples

Then you can test that the code is working by running the `quick` example.
```
# (Make sure your environment is activated)
cd examples/quick
automech run &> out.log &
```
You can see other examples in [here](./examples/).

### Tests

We cannot currently run a full test workflow on GitHub Actions because our current
set-up depends on proprietry electronic structure software.
Testing therefore involves a local workflow step, which must be done before submitting a
pull request to trigger the remaining GitHub Actions workflows.

#### Running local tests

The first time you run tests on a new machine, you will need to configure your GitHub username.
```
git config --global user.name "<username>"
```
You can then run tests locally as follows.
```
pixi run test local <node1> <node2> <...>
```
You can check the progress of this local test run as follows.
```
pixi run test status
```
You can kill the test run by simply killing all processes running on the last node,
which is where the main workflow is executed.

#### Before submitting a pull request

The above local test run is *required* before submitting a pull request with any
significant changes.
The local test data is archived and used by the GitHub Actions workflow, which also
checks the `signature.yaml` file to make sure the local tests were run with the updated
version of the code.

If you are submitting a pull request with minor changes, such as editing documentation,
that will not affect the tests, you can circumvent the local test workflow by "signing
off" on the tests with the following command.
```
pixi run test sign
```
You will be prompted to verify that the commits since the last local test workflow are
all minor and will not break the tests.
This list of untested commits will be recorded in the `signature.yaml` file, which will
override the commit hash check on GitHub Actions and allow your tests to pass.


### Subtasks

Workflow parallelization is currently not automated in AutoMech. However, if you are on a cluster with direct SSH node access and permissions to run, you can run the following commands to split an AutoMech workflow into subtasks and run them in parallel.

(1.) You can set-up these subtask jobs as follows:
```
automech subtasks setup
```
This will parse your `inp/` directory and create individual subdirectories for running each individual task for each individual species or reaction/TS. These directories will go in a folder called `subtasks/`.

(2.) If you are using the [amech-dev](https://github.com/Auto-Mech/amech-dev) Pixi environment, you can run the subtasks in parallel on a list of nodes as follows:
```
automech subtasks run csed-00{08..10}  # expands to csed-0008 csed-0009 csed-0010
```
If you are running in a different environment, you will need to pass in an activation hook as follows:
```
automech subtasks run csed-00{08..10} -a <activation hook>
```
Where the activation hook contains the bash commands to activate your environment.

(3.) To check the progress of your subtask run, you can use the following command:
```
automech subtasks status
```
This will print a color-coded table showing which tasks have failed for which species/reactions. It will also generate a `check.log` file with the paths to log files that have have not completed successfully or have a warning.

