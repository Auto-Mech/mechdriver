"""HyperQueue utilities."""

import contextlib
import math
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import TypeAlias

import jinja2
import pint
import pyparsing as pp
from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest
from hyperqueue.task.function import PythonEnv
from hyperqueue.task.task import Task
from pydantic import BaseModel, computed_field

Function: TypeAlias = Task


# Instantiate objects
def client(server_dir: str | None = None, env_prologue: str | None = None) -> Client:
    """Create HyperQueue client, for connecting to server.

    :param python_environment_prologue: Command to activate Python environment
    """
    return Client(
        server_dir=server_dir or default_server_directory(),
        python_env=PythonEnv(prologue=env_prologue or pixi_environment_prologue()),
    )


def job() -> Job:
    """Create HyperQueue job."""
    return Job()


def resource_request(cpus: int, mem: int) -> ResourceRequest:
    """Create HyperQueue resource request.

    :param cpus: Number of CPUs
    :param mem: Memory in GB
    """
    return ResourceRequest(cpus=cpus, resources={"mem": memory_mib(mem)})


# Execute system commands
def start_server(server_dir: str | None = None) -> None:
    """Start HyperQueue server."""
    # Start server and wait until it prints its info
    server_info_border_regex = re.compile(r"^\s*\+\-+\+\-+\+\s*$")
    print("Starting HyperQueue server...")
    proc = subprocess.Popen(
        ["hq", "server", "start"],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )
    # The border +---...+---...+ gets printed above and below the server info,
    # so wait until the second appearance
    print("Waiting for signal that HyperQueue server is running...")
    count = 0
    for line in proc.stdout:
        sys.stdout.write(line)
        sys.stdout.flush()
        if server_info_border_regex.match(line):
            count += 1
            if count > 1:
                break
    print("HyperQueue server started.")


def create_allocation_queue(
    mem: int, cpus: int, flags: str, manager: str, time_limit: str = "4 hr"
) -> None:
    """Create HyperQueue allocation queue."""
    time_limit_h = math.ceil(pint.Quantity(time_limit).m_as("hours"))

    # Base arguments
    args = ("hq", "alloc", "add", manager, "--time-limit", f"{time_limit_h}h")

    # Resource arguments and flags
    cpu_arg = f"--cpus={cpus}"
    mem_arg = f"--resource=mem=sum({memory_mib(mem)})"
    args += (cpu_arg, mem_arg, "--", *flags.split())

    # Extra manager-specific arguments
    if manager == "slurm":
        args += (f"--ntasks={cpus}", f"--mem-per-cpu={mem}G")

    print("HyperQueue allocation command:")
    print(" ".join(args))
    subprocess.run(args)


# Worker configuration
class WorkerConfig(BaseModel):
    """HyperQueue worker configurations."""

    name: str
    mem: int
    cpus: int = 1
    time_limit: str = "4 hr"
    idle_timeout: str = "15 min"
    manager: str | None = None
    host: str | None = None
    flags: str | None = None
    server_dir: str | None = None

    @computed_field
    @property
    def mem_mib(self) -> int:
        """Memory (MiB)."""
        return math.ceil(pint.Quantity(self.mem, "GB").m_as("MiB"))

    @computed_field
    @property
    def time_limit_hms(self) -> str:
        """Time limit (HH:MM:SS)."""
        total_s = pint.Quantity(self.time_limit).m_as("s")
        h = int(total_s // 3600)
        m = int((total_s % 3600) // 60)
        s = int(total_s % 60)
        return f"{h:02d}:{m:02d}:{s:02d}"


WORKER_DIR = Path(".workers")


def worker_configuration(
    name: str | None = None,
    mem: int | None = None,
    cpus: int | None = None,
    time_limit: str = "4 hr",
    idle_timeout: str = "15 min",
    manager: str | None = None,
    flags: str | None = None,
    host: str | None = None,
    server_dir: str | None = None,
) -> WorkerConfig:
    """Configure a worker."""
    # If `name` is None, use `host`; otherwise, set it to "worker"
    name = name or host or "worker"

    # Determine manager
    manager = determine_manager(manager)

    # Determine memory if None
    if mem is None and host is not None:
        mem = host_memory(host=host, unit="GB")

    if mem is None:
        msg = "Either memory or host must be specified."
        raise ValueError(msg)

    # Determine CPUs if None
    if cpus is None and host is not None:
        cpus = host_cpus(host=host)

    if cpus is None:
        msg = "Either CPUs or host must be specified."
        raise ValueError(msg)

    server_dir = server_dir or default_server_directory()

    # Return worker configuration
    return WorkerConfig(
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


def submit_worker_script(
    worker_config: WorkerConfig,
    template_path: str | Path | None = None,
    dry_run: bool = False,
) -> None:
    """Submit a worker script, optionally based on a template."""
    # Make sure the worker dir exists
    WORKER_DIR.mkdir(exist_ok=True)

    # Determine worker script template
    manager = worker_config.manager
    template_str = worker_script_template(manager=manager, template_path=template_path)

    # Build worker script
    template = jinja2.Template(template_str)
    script_str = template.render(**worker_config.model_dump())

    # Write worker script
    script_name = Path(worker_config.name).with_suffix(".sh")
    script_path = WORKER_DIR / script_name
    script_path.write_text(script_str)
    print(f"Script written to {script_path}")

    # Write template
    script_template_path = script_path.with_suffix(".t.sh")
    script_template_path.write_text(template_str)
    print(f"Template written to {script_template_path}")

    # Determine submission command based on manager
    if manager == "pbs":
        cmd = "qsub"
    elif manager == "slurm":
        cmd = "sbatch"
    else:
        msg = f"Workload manager '{manager}' is not a valid option ('pbs' or 'slurm')."
        raise ValueError(msg)

    # Print submission command arguments
    args = [cmd, str(script_name)]
    print(f"Submission command to be executed in {WORKER_DIR}:")
    print(" ".join(args))

    # If doing a dry run, return early
    if dry_run:
        print("Not submitting because user requested a dry run.")
        return

    # Execute submission command in worker directory
    with contextlib.chdir(WORKER_DIR):
        subprocess.run(args)


SLURM_WORKER_SCRIPT = """
#!/usr/bin/env bash
#SBATCH --nodes=1
#SBATCH --job-name={{ name }}
#SBATCH --ntasks={{ cpus }}
#SBATCH --mem={{ mem }}G
#SBATCH --time={{ time_limit_hms }}
{%- if flags is not none %}
#SBATCH {{ flags }}
{% endif %}

{%- if server_dir is not none %}
export HQ_SERVER_DIR="{{ server_dir }}"
{% endif %}

hq worker start \\
    --cpus "{{ cpus }}" \\
    --resource "mem=sum({{ mem_mib }})" \\
    --on-server-lost "stop" \\
    --idle-timeout "{{ idle_timeout }}" \\
    --time-limit "{{ time_limit }}"

"""

PBS_WORKER_SCRIPT = """
#!/usr/bin/env bash
#PBS -N {{ name }}
#PBS -l select=1:ncpus={{ cpus }}:mpiprocs={{ cpus }}
{%- if host is not none -%}
:host={{ host }}
{%- endif %}
#PBS -l walltime={{ time_limit_hms }}
{%- if flags is not none %}
#PBS {{ flags }}
{% endif %}

{%- if server_dir is not none %}
export HQ_SERVER_DIR="{{ server_dir }}"
{% endif %}

hq worker start \\
    --cpus "{{ cpus }}" \\
    --resource "mem=sum({{ mem_mib }})" \\
    --on-server-lost "stop" \\
    --idle-timeout "{{ idle_timeout }}" \\
    --time-limit "{{ time_limit }}"

"""


def worker_script_template(
    manager: str | None, template_path: str | Path | None = None
) -> str:
    """Generate a default worker script template.

    :param manager: Manually specified workload manager
    """
    if template_path is None:
        if manager == "pbs":
            return PBS_WORKER_SCRIPT
        elif manager == "slurm":
            return SLURM_WORKER_SCRIPT
        else:
            msg = f"Workload manager '{manager}' is not a valid option (pbs or slurm)"
            raise ValueError(msg)

    template_path = Path(template_path)
    return template_path.read_text()


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


def default_server_directory() -> str:
    """Path to HyperQueue server."""
    if "HQ_SERVER_DIR" in os.environ:
        return os.environ["HQ_SERVER_DIR"]

    return str(Path(os.environ["HOME"]) / ".hq-server" / "hq-current")


def pixi_environment_prologue() -> str:
    """Return Pixi Python environment prologue."""
    return subprocess.check_output(["pixi", "shell-hook"], text=True)


# Get host information over SSH
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
    mem_q = pint.Quantity(mem_expr.parse_string(mem_query.stdout).get("mem"))
    return int(mem_q.m_as(unit))


# Helpers
def memory_mib(mem: int) -> int:
    """Convert memory in GB to MiB.

    :param mem: Memory (GB)
    :return: Memory (MiB)
    """
    return math.ceil(pint.Quantity(mem, "GB").m_as("MiB"))
