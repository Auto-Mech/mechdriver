#!/usr/bin/env python
"""Local testing CLI."""

import os
import socket
import warnings

import click

import mechdriver
import test_utils as tu
from test_utils import Test


@click.group()
def main():
    """Route to test subcommand."""
    pass


@main.command("status")
def status():
    """Check the status of local tests."""
    mechdriver.subtasks.status_multiple(Test.paths())


@main.command("local")
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
def local(
    manager: str | None = None,
    time_limit: str = "4 hr",
    manager_flags: str | None = None,
    server_dir: str | None = None,
    env_prologue: str | None = None,
):
    """Run local tests on one or more nodes.

    Runs hidden local_

    :param nodes: A list of nodes
    """
    print("Process ID:", os.getpid())
    print("Host name:", socket.gethostname())

    official = False
    if os.environ.get("PIXI_ENVIRONMENT_NAME") == "test":
        official = True
        print("Official test run. Test data will be archived and committed.")
    else:
        print("Dev test run. Test data will **not** be archived or committed.")

    if manager_flags is None:
        msg = (
            "\nWARNING: Running without -f requires manual HyperQueue configuration."
            "\nMake sure you have a server running with the appropriate workers."
        )
        warnings.warn(msg, stacklevel=1)

    if manager is None:
        print("No workload manager specified with -m. Will attempt autodetection...")

    test_paths = tu.setup_tests(official=official)
    mechdriver.subtasks.setup_multiple(test_paths)
    mechdriver.subtasks.run_multiple(
        test_paths,
        manager=manager,
        time_limit=time_limit,
        manager_flags=manager_flags,
        server_dir=server_dir,
        env_prologue=env_prologue,
    )
    
    if official:
        tu.wrap_up_tests(from_archive=False, allow_override=False)


@main.command("sign")
def sign():
    """Sign off on local tests."""
    tu.wrap_up_tests(from_archive=True, allow_override=True)


if __name__ == "__main__":
    main()
