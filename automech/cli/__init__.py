import click

from . import _check_log, _run, _subtasks_run_adhoc, _subtasks_setup, _subtasks_status
from ._check_log import Status
from ._subtasks_setup import SUBTASK_DIR


@click.group()
def main():
    """AutoMech CLI"""
    pass


@main.command()
@click.option(
    "-p", "--path", default=".", show_default=True, help="The job run directory"
)
@click.option("-S", "--safemode-off", is_flag=True, help="Turn off safemode?")
def run(path: str = ".", safemode_off: bool = False):
    """Run central workflow
    Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.

    The AutoMech directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv, mechanism.dat
    """
    _run.main(path=path, safemode_off=safemode_off)


@main.command()
@click.option(
    "-p", "--path", default=".", show_default=True, help="The path to the log file"
)
def check_log(path: str = "."):
    """Check an AutoMech log file to see if it succeeded

    The path must point either directly to the log file, or to a directory where the log
    file is named "out.log"
    """
    _check_log.main(path=path)


@main.group()
def subtasks():
    """Run AutoMech subtasks in parallel"""
    pass


@subtasks.command()
@click.option(
    "-p", "--path", default=".", show_default=True, help="The job run directory"
)
@click.option(
    "-o",
    "--out-path",
    default=SUBTASK_DIR,
    show_default=True,
    help="The output path of the subtask directories",
)
@click.option(
    "-s",
    "--save-path",
    default=None,
    show_default=True,
    help="The save filesystem prefix",
)
@click.option(
    "-r",
    "--run-path",
    default=None,
    show_default=True,
    help="The run filesystem prefix",
)
@click.option(
    "-g",
    "--task-groups",
    default="els,thermo,ktp",
    show_default=True,
    help=(
        "The task groups to set up, as a comma-separated list.\n"
        "Options: els(=els-spc,els-pes), thermo, ktp"
    ),
)
def setup(
    path: str = ".",
    out_path: str = SUBTASK_DIR,
    save_path: str | None = None,
    run_path: str | None = None,
    task_groups: str = "els,thermo,ktp",
):
    """Set-up subtasks from a user-supplied AutoMech directory

    The AutoMech directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv, mechanism.dat
    """
    _subtasks_setup.main(
        path=path,
        out_path=out_path,
        save_path=save_path,
        run_path=run_path,
        task_groups=task_groups.split(","),
    )


@subtasks.command()
@click.option(
    "-p", "--path", default=SUBTASK_DIR, show_default=True, help="The subtask directory"
)
@click.option(
    "-n",
    "--nodes",
    default=None,
    show_default=True,
    help="A comma-separated list of nodes",
)
@click.option(
    "-a",
    "--activation-hook",
    default=None,
    show_default=True,
    help="An activation hook, to be called using `eval`",
)
@click.option(
    "-s",
    "--statuses",
    default=f"{Status.TBD.value}",
    show_default=True,
    help="A comma-separated list of statuses to run or re-run",
)
def run_adhoc(
    path: str = SUBTASK_DIR,
    nodes: str | None = None,
    activation_hook: str | None = None,
    statuses: str = f"{Status.TBD.value}",
):
    """Run subtasks in parallel on an Ad Hoc SSH Cluster"""
    _subtasks_run_adhoc.main(
        path=path, nodes=nodes, activation_hook=activation_hook, statuses=statuses
    )


@subtasks.command()
@click.option(
    "-p", "--path", default=SUBTASK_DIR, show_default=True, help="The subtask directory"
)
@click.option(
    "-c",
    "--check-file",
    default="check.log",
    show_default=True,
    help="Log file for writing paths to be checked",
)
@click.option(
    "-w",
    "--wrap",
    default=18,
    show_default=True,
    help="Wrap to included this many subtask columns per row",
)
def status(path: str = SUBTASK_DIR, check_file: str = "check.log", wrap: int = 18):
    """Check the status of running subtasks"""
    _subtasks_status.main(path=path, check_file=check_file, wrap=wrap)
