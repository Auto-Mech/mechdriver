#!/usr/bin/env python
"""Local testing CLI."""

import os
import socket

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
    "-f",
    "--auto-config-flags",
    default=None,
    required=True,
    help="Automatically configure HyperQueue with these sbatch/qsub flags.",
)
def local(auto_config_flags: str | None = None):
    """Run local tests on one or more nodes.

    Runs hidden local_

    :param nodes: A list of nodes
    """
    print("Process ID:", os.getpid())
    print("Host name:", socket.gethostname())

    test_paths = tu.setup_tests()
    mechdriver.subtasks.setup_multiple(test_paths)
    mechdriver.subtasks.run_multiple(test_paths, auto_config_flags=auto_config_flags)
    tu.wrap_up_tests(from_archive=False, allow_override=False)


@main.command("sign")
def sign():
    """Sign off on local tests."""
    tu.wrap_up_tests(from_archive=True, allow_override=True)


if __name__ == "__main__":
    main()
