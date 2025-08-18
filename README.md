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
3. Run the following extra command, which will install extra developer
dependencies that cannot be installed via Pixi
(see [Appendix A](#appendix-a-extra-developer-dependencies) for further
details).
```
pixi run extra-dev-installs
```
4. Check that the installation worked by running the following help command.
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
This local testing workflow is run using
[HyperQueue](https://it4innovations.github.io/hyperqueue/stable/),
which was installed with the extra developer dependencies above.
You can use the following command to confirm that the installation worked.
```
which hq
```
See [Appendix A](#appendix-a-extra-developer-dependencies) if this command was not found.

#### Running local tests

Before continuing, make sure your `git` username configuration matches your
username on GitHub
(see [Appendix B](#appendix-b-configure-your-git-username)).
You can then run tests locally as follows (see [below](#subtasks) for further details).
```
# For SLURM:
pixi run test local -m slurm -f "--partition=<partition name>"

# For PBS:
pixi run test local -m pbs -f "-q <queue name> -A <account name>"
```
The above will auto-configure a HyperQueue server to execute the testing workflow.
If the `-m` flag is ommitted, MechDriver will attempt to auto-detect Slurm or
PBS on the system.

> [!NOTE]
> On permissive clusters that allow direct SSH access to compute nodes,
> Slurm/PBS will not be aware of the resources consumed by processes executed directly over SSH.
> In this case, it may be better to manually configure the server, setting up workers to run on
> particular nodes, which can be done as follows:
> ```
> hq server start &> server.log &
> pixi run test local &> test.log &
> automech subtasks start-worker -o <node 1 name> -f "<PBS/Slurm flags>"
> automech subtasks start-worker -o <node 2 name> -f "<PBS/Slurm flags>"
> ...
> ```

You can check the progress of this local test run as follows.
```
pixi run test status
```
You can kill the test run by stopping the HyperQueue server as follows.
```
hq server stop
```

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

## Subrepos

This repository includes several submodules that also exist as separate repositories:

 - [MechAnalyzer](./src/_mechanalyzer/): Mechanism pre- and post-processing (see [here](https://github.com/Auto-Mech/mechanalyzer))
 - [AutoFile](./src/_autofile/): Filesystem databasing (see [here](https://github.com/Auto-Mech/autofile))
 - [AutoIO](./src/_autoio/): I/O interfaces to external programs (see [here](https://github.com/Auto-Mech/autoio))
 - [AutoChem](./src/_autochem/): Cheminformatics and coordinate transformation (see [here](https://github.com/Auto-Mech/autochem))

If you wish to be able to sync these submodules with their remote repositories,
this can be done using
[Git Subrepo](https://github.com/ingydotnet/git-subrepo?tab=readme-ov-file#commands),
which was installed with the extra developer dependencies above.
You can use the following command to confirm that the installation worked.
```
git subrepo -h
```
See [Appendix A](#appendix-a-extra-developer-dependencies) if this command was not found.

### Syncing

Before continuing, make sure your `git` username configuration matches your
username on GitHub
(see [Appendix B](#appendix-b-configure-your-git-username)).

To pull updates for one or more subrepos, you can use the `pull` task.
```
pixi run pull all     # pull changes for all subrepos
pixi run pull autoio  # pull changes for AutoIO only
```
To push updates back to the subrepos, you can use the `push` task.
```
pixi run push all     # push changes for all subrepos
pixi run push autoio  # push changes for AutoIO only
```

### Advanced

The above Pixi tasks are sufficient for working with the default branches of each fork
and keeping them in sync with their upstream repositories.
To pull from/push to a specific branch of a subrepo, you can add a `-b` flag.
```
pixi run pull autoio -b <branch name>
pixi run push autoio -b <branch name>
```
This flag, along with any others added after the repository name, is simply passed along
to the `git subrepo pull` and `git subrepo push` commands, which are documented
[here](https://github.com/ingydotnet/git-subrepo?tab=readme-ov-file#commands).

For more advanced usage, you can run the commands manually, which is facilitated by the
tab completion provided by `git-subrepo`.
For example, the pull/push commands above are equivalent to the following.
```
git subrepo pull src/_autoio -r git@github.com:<username>/autoio.git -b <branch name>
git subrepo pull src/_autoio -r git@github.com:<username>/autoio.git -b <branch name>
```
Here, the [GitHub CLI](https://cli.github.com/manual/) can also come in handy. It is
installed in the `dev` Pixi environment of this repository.
For example, you can sync a particular branch of your fork against its uptream [as follows](https://cli.github.com/manual/gh_repo_sync).
```
gh repo sync <username>/autoio -b <branch name>
```

## Subtasks

Workflow parallelization is currently not fully automated in AutoMech. However, you can
run the following commands to split an AutoMech workflow into subtasks and then run them
using [HyperQueue](https://it4innovations.github.io/hyperqueue/stable/).
See [Appendix A](#appendix-a-extra-developer-dependencies) for instructions on
installing HyperQueue.

To see if it works, you can test the following steps on the same
["quick" example](examples/quick/) that you ran above.

**Setup.** You can set up these subtask jobs as follows:
```
automech subtasks setup
```
This will parse your `inp/` directory and create individual subdirectories for running each individual task for each individual species or reaction/TS. These directories will go in a folder called `subtasks/`.

**Run.** You can run these subtask jobs as follows.
```
# For SLURM:
automech subtasks run -f "--partition=<partition name>"

# For PBS:
automech subtasks run -f "-q <queue name> -A <account name>"
```
The `-f` flag allows you to pass additional flags to SLURM or PBS. This would be
anything beyond basic resources (memory, CPUs, etc.) that you are required to put in
your `sbatch` or `qsub` scripts.

**Check status.** To check the progress of your subtask run, you can use the following command:
```
automech subtasks status
```
This will print a color-coded table showing which tasks have failed for which species/reactions. It will also generate a `check.log` file with the paths to log files that have have not completed successfully or have a warning.


## Appendix A: Extra Developer Dependencies

The simplest way to install the extra developer dependencies is by running the
following Pixi task.
```
pixi run extra-dev-installs
```
If you run into issues with this, you can manually install the two dependencies,
HyperQueue and Git Subrepo, as follows.

**HyperQueue**
HyperQueue allows you to execute parallel workflows on PBS or SLURM.
You can install it as follows.
```
wget https://github.com/It4innovations/hyperqueue/releases/download/v0.20.0/hq-v0.20.0-linux-x64.tar.gz
tar -zxvf hq-v0.20.0-linux-x64.tar.gz -C /directory/in/shell/path
```
The second command puts the HyperQueue executable into a directory that is in your shell
path.  For example, this might be `$HOME/bin` if you have `export PATH=$PATH:$HOME/bin`
in your `.bashrc`.


**Git Subrepo**
Git Subrepo allows you to sync the submodules in this repository with their
remote counterparts.
You can install it as follows.
```
git clone https://github.com/ingydotnet/git-subrepo /path/to/git-subrepo
echo 'source /path/to/git-subrepo/.rc' >> ~/.bashrc
```


## Appendix B: Configure Your Git Username

The first time you run tests on a new machine, you will need to configure your
username with `git`.
```
git config --global user.name "<username>"
```
Make sure this matches your username on GitHub.
