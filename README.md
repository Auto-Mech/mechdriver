# MechDriver

MechDriver is the main driver for executing AutoMech workflows.

## Install

MechDriver can be installed by your favorite conda package manager from the
[`auto-mech` channel](https://anaconda.org/Auto-Mech/mechdriver), with
dependencies from the `conda-forge` channel.
We would highly recommend the Pixi package manager, which you can easily install
with [this command](https://pixi.sh/latest/#installation).
The Pixi command to install MechDriver is as follows.
```
pixi global install mechdriver -c conda-forge -c auto-mech
```
This "global install" command makes the code available across all environments.
Like Pixi itself, however, it only uses your home directory and does not require
special permissions, so you can do this anywhere.
If you prefer, Pixi also allows you to install the code into
[a local workspace environment](https://pixi.sh/latest/first_workspace/).

If the installation worked, you should have `mechdriver` and `mechanalyzer` in
your path.
You can check this by running the following help commands.
```
mechdriver --help
mechanalyzer --help
```
Note that, if you are not using a globally installed version of the code, you
will need to activate the appropriate environment first.

## Run

The easiest way to start running MechDriver is to try out one of the test
examples
These can be downloaded as follows:
```
wget https://github.com/avcopan/mechdriver/raw/refs/heads/dev/tests/archive.tgz
mkdir examples
tar -xzf archive.tgz -C examples
```

You can then run the `quick` example as follows.
```
cd examples/quick
mechdriver run
```
This will run very quickly because it comes with pre-populated electronic
structure data.

> [!TIP]
> To see live output as the calculation runs, you may need to turn off Python
> output buffering:
> ```
> export PYTHONUNBUFFERED=1
> ```
> You can put this line in your `.bashrc`.

To run from scratch, make sure the `g16` executable for Gaussian 16 is in your path.
Then remove the `save` directory for the example and try running it again.
This time, you may wish to write the output to a log file.
```
mechdriver run &> out.log &
```
Electronic structure calculations will run in a directory called `run`.
The paths to these calculations will be printed to the log file as the workflow
proceeds.

## Developer Instructions

### Setup

> [!IMPORTANT]
> Until you complete the developer set-up, you will need to add the `--frozen`
> flag to all Pixi commands.

To get started, you will need to have
[installed Pixi](https://pixi.sh/latest/#installation)
and forked the repositories as described in
[Appendix A](#appendix-a-forking-submodule-repositories).
You will also need to have configured your Git username as described in
[Appendix B](#appendix-b-configure-your-git-username).

You can then complete your developer set-up by cloning this repository and
running the following command inside it.
```
pixi run --frozen dev-setup
```
Follow the prompts to configure the set-up to your liking.
Among other configurations, this set-up clones the other Auto-Mech repositories
into the parent directory of your mechdriver repository.
```
ls ..
autochem  autofile  autoio  mechanalyzer  mechdriver
```
If you run into issues, see [Appendix C](#appendix-c-manual-developer-setup) for
manual set-up instructions.

After re-starting your shell, you can test that the setup worked by activating
the development environment and running the help command.
```
mechenv
mechdriver --help
```
The `mechenv` command is a bash function for activating MechDriver environments
that was added to your `.bashrc` configuration by the developer set-up.
It activates the `dev` environment by default, but you can also pass it the name
of a another environment, such as `mechenv test` to activate the `test`
environment.
For more information on activating Pixi environments, see
[here](https://pixi.sh/dev/workspace/environment/#activation).

You can then run any of the examples [in the `examples/` directory](./examples/)
as [described above](#run).

### Git Helper Tasks

The following Pixi tasks are available to facilitate Git operations with your
five local Auto-Mech repositories:

1. `pixi run pull`. This does a `git pull --rebase upstream dev` on all five
repositories to update them against the central Auto-Mech upstream.
2. `pixi run git <commands>`. This runs whatever `git` commands you pass it in
all five repositories.

For example, a common operation you might wish to do is:
```
pixi run pull                   # Rebase each repo against upstream
pixi run git push origin dev    # Push each repo to origin
```

### Test

> [!TIP]
> This testing workflow uses the experimental "subtasks" module in MechDriver to
> parallelize subtask execution with
> [HyperQueue](https://it4innovations.github.io/hyperqueue/stable/).
> See [below](#subtasks) for more information on subtask parallelization.

The main MechDriver test workflow must be run locally on a cluster with access
to Gaussian 16.
The archived data from this local test run are then used for automated testing on
GitHub Actions.

Assuming you have completed the developer set-up, you can initiate the test
workflow as follows.
```
# Activate test environment
mechenv test

# Run test workflow via SLURM ...
pixi run test local -m slurm -f "--partition=<partition name>"

# ... OR run test workflow via PBS
pixi run test local -m pbs -f "-q <queue name> -A <account name>"
```
The above will auto-configure HyperQueue to execute the testing workflow.
In some cases, you may wish to configure HyperQueue to run on particular nodes.
This can be done as described in
[Appendix D](#appendix-d-running-test-workflow-on-select-nodes).

> [!IMPORTANT]
> The `test` environment uses the **published conda packages** for `autochem`,
> `autoio`, `autofile`, and `mechanalyzer` specified in `pyproject.toml`.
> These can be updated
> [as described below](#updating-conda-packages).
> To test against your local changes to these repositories before updating the
> conda packages, you can activate the `dev` environment.
> In this case, the archive for GitHub Actions will not be generated.

You can check the progress of a local test workflow as follows.
```
pixi run test status
```
You can kill a local test workflow by stopping the HyperQueue server as follows.
```
hq server stop
```

### Signing Tests

The above local test run is *required* before submitting a pull request with any
significant changes.
The local test data is archived and used by the GitHub Actions workflow, which also
checks the `signature.yaml` file to make sure the local tests were run against
the incoming version of the code.

If you are submitting a pull request with minor changes, such as editing documentation,
that will not affect the tests, you can circumvent the local test workflow by "signing
off" on the tests with the following command.
```
pixi run test sign
```
You will be prompted to affirm that the commits since the last local test
workflow are all minor and will not break the tests.
This list of untested commits will be recorded in the `signature.yaml` file, which will
override the commit hash check on GitHub Actions and allow your tests to pass.

### Updating Conda Packages

> [!WARNING]
> Versions are now automatically handled by the BumpVer versioning tool.
> Therefore, **do not** manually change the version of a given package in its
> `pyproject.toml`.  Instead, the version will be automatically updated by
> triggering the release workflow as described below.

To update the conda package for one of the AutoMech repositories, you can
trigger its `release` workflow by following these steps:

1. Submit a pull request as usual.
2. In the pull request UI, click "Labels" and add the label `release:patch`.
3. If and **only if** the tests pass, merge the pull request.

The `release:patch` label is the signal that tells GitHub Actions to run the
`release` workflow.
This will (1.) bump the version number, (2.) publish the conda package to the
[`auto-mech` channel](https://anaconda.org/Auto-Mech/), and (3.) create an
associated GitHub release.

If you have done this for one or more of the lower-level repositories, you will
also need to update the MechDriver `pyproject.toml` file with their new version
numbers.
You can do so as follows:
```
pixi run pull   # Pull the release commits from upstream
pixi run update # Update the version numbers in pyproject.toml
```
This will update the `package.run-dependencies` and `dependencies` tables in
your `pyproject.toml` with the new version numbers of the lower-level
repositories and update the lockfile.
If you make any manual changes to these tables,
**make sure that they match exactly**.
This is necessary to ensure that the packaged version of the code matches the
tested version.

## Experimental Feature: Subtask Parallelization

As an experimental feature in MechDriver, you can parallelize across subtasks in
your workflow with
[HyperQueue](https://it4innovations.github.io/hyperqueue/stable/).
If you are not using a developer set-up, you will first need to
[install the HyperQueue executable](https://it4innovations.github.io/hyperqueue/stable/installation/#binary-distribution-recommended)
and
[the HyperQueue Python API](https://it4innovations.github.io/hyperqueue/stable/python/).

To try it out, you can test the following steps on the same
["quick" example](examples/quick/) that you ran above.

**Setup.** You can set up these subtask jobs as follows:
```
mechdriver subtasks setup
```
This will parse your `inp/` directory and create individual subdirectories for running each individual task for each individual species or reaction. These directories will go in a folder called `subtasks/`.

**Run.** You can run these subtask jobs as follows.
```
# For SLURM:
mechdriver subtasks run -f "--partition=<partition name>"

# For PBS:
mechdriver subtasks run -f "-q <queue name> -A <account name>"
```
The `-f` flag allows you to pass additional flags to SLURM or PBS. This would be
anything beyond basic resources (memory, CPUs, etc.) that you are required to put in
your `sbatch` or `qsub` scripts.

**Check status.** To check the progress of your subtask run, you can use the following command:
```
mechdriver subtasks status
```
This will print a color-coded table showing which tasks have failed for which species/reactions. It will also generate a `check.log` file with the paths to log files that have have not completed successfully or have a warning.


## Appendix A: Forking the AutoMech Repositories

To get started as a new developer, log into GitHub (or create an account) and follow
[these instructions](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo#forking-a-repository)
to fork the following five repositories:
 - [AutoChem](https://github.com/Auto-Mech/autochem)
 - [AutoIO](https://github.com/Auto-Mech/autoio)
 - [AutoFile](https://github.com/Auto-Mech/autofile)
 - [MechAnalyzer](https://github.com/Auto-Mech/mechanalyzer)
 - [MechDriver](https://github.com/Auto-Mech/mechdriver)

## Appendix B: Configuring Your Git Username

The first time getting set up on a new machine, you will need to configure your
username with `git`.
```
git config --global user.name "<username>"
```
Make sure this matches your username on GitHub.

## Appendix C: Manual Developer Setup

If you run into issues with the automated developer setup described above, you
can get set up manually as follows.

1. On GitHub, follow
[these instructions](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork#syncing-a-fork-branch-from-the-web-ui)
to make sure all of your forks are synced.
Then clone your fork for each of repository into the parent directory containing
this `mechdriver/` repository.
```
cd ..
git clone git@github.com:<username>/autochem.git
git clone git@github.com:<username>/autoio.git
git clone git@github.com:<username>/autofile.git
git clone git@github.com:<username>/mechanalyzer.git
```
2. Download the HyperQueue executable to a directory of your choosing.
```
wget https://github.com/It4innovations/hyperqueue/releases/download/v0.20.0/hq-v0.20.0-linux-x64.tar.gz
tar -zxvf hq-v0.20.0-linux-x64.tar.gz -C $HOME/bin
```
3. In your `.bashrc` file, add the directory with the HyperQueue executable to
your path, set an environment variable to turn off Python buffering, and add a
bash function for easily activating this Pixi environment from anywhere on your
system.
```
# In .bashrc:
export PATH="$HOME/bin:PATH"
export PYTHONUNBUFFERED=1
mechenv() {
    local env_name="${1:-dev}"
    eval "$(pixi shell-hook -e "$env_name" --manifest-path "/path/to/mechdriver")"
}
```

## Appendix D: Running Test Workflow on Select Nodes

On clusters where users frequently run without the Slurm or PBS workload
manager, you may need to instead configure HyperQueue to run the test workflow on particular nodes.
This can be done as follows.
```
hq server start &> server.log &
pixi run test local &> test.log &
mechdriver subtasks start-worker -o <node 1 name> -f "<PBS/Slurm flags>"
mechdriver subtasks start-worker -o <node 2 name> -f "<PBS/Slurm flags>"
...
```
