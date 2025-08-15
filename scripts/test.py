#!/usr/bin/env python
"""Local testing CLI."""

import contextlib
import os
import shutil
import socket
import subprocess
import warnings
from pathlib import Path

import click
import pint
import pyparsing as pp

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
    "-f",
    "--manager-flags",
    default=None,
    help=(
        "Automatically configure HyperQueue allocation with "
        "these PBS or Slurm submission flags"
    ),
)
def local(manager: str | None = None, manager_flags: str | None = None):
    """Run local tests on one or more nodes.

    Runs hidden local_

    :param nodes: A list of nodes
    """
    print("Process ID:", os.getpid())
    print("Host name:", socket.gethostname())

    if manager_flags is None:
        msg = (
            "\nWARNING: Running without -f requires manual HyperQueue configuration."
            "\nMake sure you have a server running with the appropriate workers."
        )
        warnings.warn(msg, stacklevel=1)

    if manager is None:
        print("No workload manager specified with -m. Will attempt autodetection...")

    test_paths = tu.setup_tests()
    mechdriver.subtasks.setup_multiple(test_paths)
    mechdriver.subtasks.run_multiple(
        test_paths,
        manager=manager,
        manager_flags=manager_flags,
    )
    tu.wrap_up_tests(from_archive=False, allow_override=False)


@main.command("sign")
def sign():
    """Sign off on local tests."""
    tu.wrap_up_tests(from_archive=True, allow_override=True)


# Create node worker
WORKER_DIR = Path(".workers")


@main.command("create-node-worker")
@click.argument("host")
@click.option("-q", "--queue", required=True, help="Queue name")
@click.option("-A", "--account", required=True, help="Account name")
def create_node_worker(host: str, queue: str, account: str):
    """Create a worker for a specific node."""
    # Create script
    cpus = host_cpus(host)
    mem_mib = host_memory(host, unit="MiB")
    name = f"worker-{host}"
    script_text = worker_script_text(
        host, name=name, cpus=cpus, mem_mib=mem_mib, queue=queue, account=account
    )

    # Submit script
    WORKER_DIR.mkdir(exist_ok=True)
    with contextlib.chdir(WORKER_DIR):
        script_name = f"{name}.sh"
        script = Path(script_name)
        script.write_text(script_text)
        subprocess.run(["qsub", script_name])


# Helper functions
def host_cpus(host: str) -> int:
    """Determine host number of CPUs."""
    cpus_query = subprocess.run(
        ["ssh", host, "nproc --all"], capture_output=True, text=True
    )
    return int(cpus_query.stdout)


def host_memory(host: str, unit: str = "GB") -> int:
    """Determine host memory."""
    mem_query = subprocess.run(
        ["ssh", host, "grep MemTotal /proc/meminfo"], capture_output=True, text=True
    )
    mem_expr = pp.Suppress("MemTotal:") + pp.SkipTo(pp.StringEnd())("mem")
    mem = pint.Quantity(mem_expr.parse_string(mem_query.stdout).get("mem"))
    return int(mem.m_as(unit))


PBS_SCRIPT = """
#!/bin/bash
#PBS -l select=1:host={host}
#PBS -N {name}
#PBS -l walltime=01:00:00
#PBS -q {queue} -A {account}

hq worker start \\
    --idle-timeout "1h" \\
    --manager "pbs" \\
    --cpus "{cpus:d}" \\
    --resource "mem=sum({mem_mib})" \\
    --on-server-lost "stop" \\
    --time-limit "2h"

"""


def worker_script_text(
    host: str, *, name: str, cpus: int, mem_mib: int, queue: str, account: str
) -> str:
    """Generate worker script text."""
    if shutil.which("qsub"):
        return PBS_SCRIPT.format(
            host=host,
            name=name,
            queue=queue,
            account=account,
            cpus=cpus,
            mem_mib=mem_mib,
        )

    msg = "Worker script generation currently only implemented for PBS."
    raise NotImplementedError(msg)


if __name__ == "__main__":
    main()
