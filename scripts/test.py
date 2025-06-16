#!/usr/bin/env python
"""Local testing CLI."""
import os
import socket
import subprocess
from collections.abc import Sequence

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


@main.command("local", hidden=True)
@click.argument("nodes", nargs=-1)
def local(nodes: Sequence[str]):
    """Run local tests on one or more nodes.

    Runs hidden local_

    :param nodes: A list of nodes
    """
    log_name = "test.log"
    test_cmd = " ".join(["pixi run test local_", *nodes])
    cmd = ["pixi", "run", "node", nodes[-1], log_name, test_cmd]
    print(subprocess.check_output(cmd, text=True))


@main.command("local_", hidden=True)
@click.argument("nodes", nargs=-1)
def local_(nodes: Sequence[str]):
    """Run local tests on one or more nodes.

    :param nodes: A list of nodes
    """
    print("Process ID:", os.getpid())
    print("Host name:", socket.gethostname())

    test_paths = tu.setup_tests()
    mechdriver.subtasks.setup_multiple(test_paths)
    mechdriver.subtasks.run_multiple(
        test_paths, nodes=nodes, activation_hook=tu.pixi_activation_hook()
    )
    tu.wrap_up_tests(from_archive=False, allow_override=False)


@main.command("sign")
def sign():
    """Sign off on local tests."""
    tu.wrap_up_tests(from_archive=True, allow_override=True)


if __name__ == "__main__":
    main()
