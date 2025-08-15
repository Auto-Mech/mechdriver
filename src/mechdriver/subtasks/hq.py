"""HyperQueue utilities."""

import math
import os
import shutil
import subprocess
import time
from pathlib import Path
from typing import TypeAlias

import pint
from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest
from hyperqueue.task.function import PythonEnv
from hyperqueue.task.task import Task

Function: TypeAlias = Task


# Instantiate objects
def client(env_prologue: str | None = None) -> Client:
    """Create HyperQueue client, for connecting to server.

    :param python_environment_prologue: Command to activate Python environment
    """
    return Client(
        server_dir=current_server_path(),
        python_env=PythonEnv(prologue=env_prologue or pixi_environment_prologue()),
    )


def job() -> Job:
    """Create HyperQueue job."""
    return Job()


def resource_request(cpus: int, mem: int) -> ResourceRequest:
    """Create HyperQueue resource request."""
    return ResourceRequest(cpus=cpus, resources={"mem": memory_mib(mem)})


# Execute system commands
def start_server() -> None:
    """Start HyperQueue server."""
    server_path = current_server_path()
    subprocess.Popen(["hq", "server", "start"])
    # Wait up to 1 second for the file to appear
    for _ in range(10):
        time.sleep(0.1)
        if os.path.exists(server_path):
            break
    assert os.path.exists(server_path), f"Could not start server at {server_path}"


def create_allocation_queue(
    mem: int, cpus: int, flags: str, manager: str | None = None
) -> None:
    """Create HyperQueue allocation queue."""
    # Determine workload manager
    manager = determine_manager(manager=manager)

    # Base arguments
    args = ("hq", "alloc", "add", manager, "--time-limit", "4h")

    # Resource arguments and flags
    cpu_arg = f"--cpus={cpus}"
    mem_arg = f"--resource=mem=sum({memory_mib(mem)})"
    args += (cpu_arg, mem_arg, "--", *flags.split())

    # Extra manager-specific arguments
    if manager == "slurm":
        args += ("--ntasks=1", f"--mem={mem}G")

    print("HyperQueue allocation command:")
    print(" ".join(args))
    subprocess.run(args)


# Get system information
def determine_manager(manager: str | None = None) -> str:
    """Detect which workload manager (PBS or Slurm) is on the system.

    :param manager: Manually specified workload manager
    """
    if manager is None:
        if shutil.which("qsub"):
            manager = "pbs"
        elif shutil.which("sbatch"):
            manager = "slurm"

    if manager is None:
        msg = "No SLURM or PBS detected. Please manually configure HyperQueue."
        raise ValueError(msg)

    manager = manager.lower()

    if manager not in ("pbs", "slurm"):
        msg = f"Workload manager '{manager}' is not a valid option ('pbs' or 'slurm')."
        raise ValueError(msg)

    return manager


def current_server_path() -> Path:
    """Path to HyperQueue server."""
    return Path(os.environ["HOME"]) / ".hq-server" / "hq-current"


def pixi_environment_prologue() -> str:
    """Return Pixi Python environment prologue."""
    return subprocess.check_output(["pixi", "shell-hook"], text=True)


# Helpers
def memory_mib(mem: int) -> int:
    """Convert memory in GB to MiB.

    :param mem: Memory (GB)
    :return: Memory (MiB)
    """
    return math.ceil(pint.Quantity(mem, "GB").m_as("MiB"))
