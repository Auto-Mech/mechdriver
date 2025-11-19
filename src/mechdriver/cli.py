import warnings
from collections.abc import Sequence
from pathlib import Path

import click

from . import subtasks, tools
from .base import Status, check_log, run
from .subtasks import hq


@click.group()
def main():
    """MechDriver CLI"""
    pass


@main.command("run")
@click.option(
    "-p", "--path", default=".", show_default=True, help="The job run directory"
)
@click.option("-S", "--safemode-off", is_flag=True, help="Turn off safemode?")
def run_(path: str = ".", safemode_off: bool = False):
    """Run central workflow.

    Central Execution script to launch a MechDriver process which will
    parse all of the user-supplied input files in a specified directory, then
    launches all of the requested electronic structure, transport,
    thermochemistry and kinetics calculations via their associated
    sub-drivers.

    The AutoMech directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv, mechanism.dat
    """
    run(path=path, safemode_off=safemode_off)


@main.command("check-log")
@click.option(
    "-p", "--path", default=".", show_default=True, help="The path to the log file"
)
def check_log_(path: str = "."):
    """Check an AutoMech log file to see if it succeeded

    The path must point either directly to the log file, or to a directory where the log
    file is named "out.log"
    """
    check_log(path=path, log=True)


@main.group("subtasks")
def subtasks_():
    """Run AutoMech subtasks in parallel"""
    pass


@subtasks_.command("setup")
@click.argument("paths", nargs=-1)
@click.option(
    "-n",
    "--dir-name",
    default=subtasks.SUBTASK_DIR,
    show_default=True,
    help="The subtask directory name",
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
    default="els,thermo,ktp,proc",
    show_default=True,
    help=(
        "The task groups to set up, as a comma-separated list.\n"
        "Options: els(=els-spc,els-pes), thermo, ktp"
    ),
)
def subtasks_setup_(
    paths: tuple[str, ...],
    dir_name: str = subtasks.SUBTASK_DIR,
    save_path: str | None = None,
    run_path: str | None = None,
    task_groups: str = "els,thermo,ktp",
):
    """Set up subtasks from a user-supplied AutoMech directory

    Passing multiple paths will do this for multiple directories.

    Each AutoMech directory must contain an `inp/` subdirectory with the following
    required files: run.dat, theory.dat, models.dat, species.csv, mechanism.dat
    """
    paths = paths if paths else (".",)
    subtasks.setup_multiple(
        paths=paths,
        dir_name=dir_name,
        save_path=save_path,
        run_path=run_path,
        task_group_keys=task_groups.split(","),
    )


@subtasks_.command("run")
@click.argument("paths", nargs=-1)
@click.option(
    "-d",
    "--dir-name",
    default=subtasks.SUBTASK_DIR,
    show_default=True,
    help="The subtask directory name",
)
@click.option(
    "-t",
    "--statuses",
    default=f"{Status.TBD.value}",
    show_default=True,
    help="A comma-separated list of statuses to run or re-run",
)
@click.option(
    "-m",
    "--manager",
    default=None,
    help="Specify workload manager (PBS or Slurm) instead of autodetecting",
)
@click.option(
    "-l",
    "--time-limit",
    default="4 hr",
    show_default=True,
    help="Worker time limit with units, e.g. '4 hr'",
)
@click.option(
    "-f",
    "--manager-flags",
    default=None,
    help=(
        "Automatically configure HyperQueue allocation with "
        "these PBS or Slurm submission flags"
    ),
)
@click.option(
    "-s",
    "--server-dir",
    default=None,
    help="HyperQueue server directory",
)
@click.option(
    "-e",
    "--env-prologue",
    default=None,
    help="Command(s) to activate Python environment",
)
def subtasks_run_(
    paths: Sequence[str] = (".",),
    dir_name: str = subtasks.SUBTASK_DIR,
    statuses: str = f"{Status.TBD.value}",
    time_limit: str = "4 hr",
    manager: str | None = None,
    manager_flags: str | None = None,
    server_dir: str | None = None,
    env_prologue: str | None = None,
):
    """Run subtasks in parallel using HyperQueue.

    Use a space-separated list of paths
        path1 path2 path3
    or
        path{1..3}

    """
    if manager_flags is None:
        msg = (
            "\nWARNING: Running without -f requires manual HyperQueue configuration."
            "\nMake sure you have a server running with the appropriate workers."
        )
        warnings.warn(msg, stacklevel=1)

    if manager is None:
        print("No workload manager specified with -m. Will attempt autodetection...")

    paths = paths if paths else (".",)
    subtasks.run_multiple(
        paths=paths,
        dir_name=dir_name,
        statuses=list(map(Status, statuses.split(","))),
        time_limit=time_limit,
        manager=manager,
        manager_flags=manager_flags,
        server_dir=server_dir,
        env_prologue=env_prologue,
    )


@subtasks_.command("status")
@click.argument("paths", nargs=-1)
@click.option(
    "-n",
    "--dir-name",
    default=subtasks.SUBTASK_DIR,
    show_default=True,
    help="The subtask directory name",
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
def subtasks_status_(
    paths: tuple[str, ...],
    dir_name=subtasks.SUBTASK_DIR,
    check_file: str = "check.log",
    wrap: int = 18,
):
    """Check the status of running subtasks"""
    paths = paths if paths else (".",)
    subtasks.status_multiple(
        paths=paths, dir_name=dir_name, check_file=check_file, wrap=wrap
    )


@subtasks_.command("start-worker")
@click.option(
    "-t",
    "--template",
    default=None,
    help="Script template",
)
@click.option(
    "-n",
    "--name",
    default=None,
    help="Worker name",
)
@click.option(
    "-m",
    "--mem",
    default=None,
    show_default=True,
    help="Memory in GB.",
)
@click.option(
    "-c",
    "--cpus",
    default=None,
    show_default=True,
    help="Number of CPUs",
)
@click.option(
    "-l",
    "--time-limit",
    default="4 hr",
    show_default=True,
    help="Worker time limit with units, e.g. '4 hr'",
)
@click.option(
    "-i",
    "--idle-timeout",
    default="15 min",
    show_default=True,
    help="Worker idle timeout with units, e.g. '15 min'",
)
@click.option(
    "-w",
    "--manager",
    default=None,
    help="Specify workload manager (PBS or Slurm) instead of autodetecting",
)
@click.option(
    "-f",
    "--flags",
    default=None,
    help="PBS or Slurm submission flags",
)
@click.option(
    "-o",
    "--host",
    default=None,
    help="Specify a host name",
)
@click.option(
    "-s",
    "--server-dir",
    default=None,
    help="HyperQueue server directory",
)
@click.option(
    "-d", "--dry-run", is_flag=True, help="Dry run: Generate script without submitting."
)
def start_worker(
    template: str | None = None,
    name: str | None = None,
    mem: int | None = None,
    cpus: int | None = None,
    time_limit: str = "4 hr",
    idle_timeout: str = "15 min",
    manager: str | None = None,
    flags: str | None = None,
    host: str | None = None,
    server_dir: str | None = None,
    dry_run: bool = False,
):
    """Start a HyperQueue worker for running subtasks."""
    worker_config = hq.worker_configuration(
        name=name,
        cpus=cpus,
        mem=mem,
        time_limit=time_limit,
        idle_timeout=idle_timeout,
        manager=manager,
        flags=flags,
        host=host,
        server_dir=server_dir,
    )
    hq.submit_worker_script(worker_config, template_path=template, dry_run=dry_run)


@main.command()
@click.option(
    "-o",
    "--options",
    default="insert_options.txt",
    show_default=True,
    help="option file, example in mechdriver/automech/tools/insert_options.txt",
)
def automated_insert(options):
    insert_options_dct = tools.automated_insert.parse_script_options(options)
    tools.automated_insert.main(insert_options_dct)
