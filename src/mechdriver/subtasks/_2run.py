"""Standalone script to run AutoMech subtasks in parallel using HyperQueue."""

import contextlib
import functools
import itertools
import math
import os
import shutil
import subprocess
import time
from collections.abc import Callable, Sequence
from pathlib import Path

import networkx as nx
import pint
import yaml
from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest
from hyperqueue.task.function import PythonEnv
from hyperqueue.task.task import Task as HQTask

from ..base import Extension, Status
from ..base import run as run_automech
from ._0setup import INFO_FILE, SUBTASK_DIR, SubtasksInfo, Task
from ._1status import log_paths_with_check_results, parse_subtask_status

HQTaskKey = tuple[int, int, str]

HOME = Path(os.environ["HOME"])
HQ_PATH = HOME / ".hq-server" / "hq-current"

SCRIPT_DIR = Path(__file__).parent / "scripts"


class Script:
    ignore_error = str(SCRIPT_DIR / "ignore_error.sh")
    lock_file = str(SCRIPT_DIR / "lock_file.sh")


def run_multiple(
    paths: Sequence[str | Path] = (".",),
    dir_name: str = SUBTASK_DIR,
    statuses: Sequence[Status] = (Status.TBD,),
    auto_config_flags: str | None = None,
    python_environment: str | None = None,
) -> None:
    """Run multiple sets of subtasks in parallel using HyperQueue.

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param paths: The paths where the AutoMech subtasks were set up
    :param dir_name: The subtask directory name
    :param hyperqueue_path: The path to the HyperQueue server directory
    :param statuses: A comma-separated list of status to run or re-run
    :param auto_config_flags: Sbatch/qsub flags for HyperQueue autoconfiguration
    :param python_environment: Command to activate Python environment
    """
    if auto_config_flags is not None:
        start_hyperqueue_server()

    # Set up the HyperQueue client
    python_environment = python_environment or subprocess.check_output(
        ["pixi", "shell-hook"], text=True
    )

    client = Client(HQ_PATH, python_env=PythonEnv(prologue=python_environment))

    # Set up the HyperQueue job workflow
    job = Job()

    for path in paths:
        job = setup_job(
            path=path,
            dir_name=dir_name,
            statuses=statuses,
            auto_config_flags=auto_config_flags,
            job=job,
        )

    submitted_job = client.submit(job)
    client.wait_for_jobs([submitted_job])


def setup_job(
    path: str | Path = ".",
    dir_name: str = SUBTASK_DIR,
    statuses: Sequence[Status] = (Status.TBD,),
    auto_config_flags: str | None = None,
    job: Job | None = None,
) -> Job:
    """Run subtasks in parallel using HyperQueue.

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param dir_name: The subtask directory name
    :param statuses: A comma-separated list of status to run or re-run
    :param auto_config: Automatically configure HyperQueue with these sbatch/qsub flags
    :param job: Append to an existing job
    """
    path = Path(path).resolve()
    dir_path = path / dir_name
    info_file = dir_path / INFO_FILE
    info = SubtasksInfo.model_validate(yaml.safe_load(info_file.read_text()))

    if auto_config_flags is not None:
        all_tasks = list(itertools.chain.from_iterable(info.task_groups))
        mem = max(t.mem for t in all_tasks)
        cpus = max(t.nprocs for t in all_tasks)
        add_hyperqueue_allocation(mem=mem, cpus=cpus, flags=auto_config_flags)

    # Make sure the run and save directories exist
    info.run_path.mkdir(exist_ok=True)
    info.save_path.mkdir(exist_ok=True)

    # Set up the HyperQueue job workflow
    job = job or Job()

    # For now, just do this for the first task group
    hq_task_dct: dict[HQTaskKey, HQTask] = {}
    dep_graph = dependency_graph(info.task_groups)
    for group_idx, task_group in enumerate(info.task_groups):
        for task_idx, task in enumerate(task_group):
            for subtask in task.subtasks:
                # Determine dependencies from dependency graph
                hq_task_key = (group_idx, task_idx, subtask.key)
                dep_hq_task_keys = dep_graph.predecessors(hq_task_key)
                dep_hq_tasks = [hq_task_dct[k] for k in dep_hq_task_keys]

                hq_task = automech_hyperqueue_task(
                    job=job,
                    path=dir_path / subtask.path,
                    log_path=dir_path / subtask.path / "out.log",
                    deps=dep_hq_tasks,
                    cpus=task.nprocs,
                    mem=task.mem,
                    workers=subtask.nworkers,
                )

                # Add the job to the job dictionary
                hq_task_dct[hq_task_key] = hq_task

    return job


def automech_hyperqueue_task(
    job: Job,
    path: Path,
    log_path: Path,
    deps: Sequence[HQTask],
    cpus: int,
    mem: int,
    workers: int = 1,
) -> HQTask:
    r"""Create a HyperQueue task to run automech.

    When n > 1 workers are requested, this creates n instances of the task with the
    original dependencies, ignoring any errors that occur. It then creates a final
    instance of the task with these n tasks as dependencies, which does *not* ignore
    errors. This serves to determine whether the n-worker task succeeded.

    So the resulting task subgraph looks like this:

                previous dependencies
               /        |      ...   \
            worker 1  worker 2 ...  worker n   <= ignore errors
               \        |      ...   /
             final task to confirm success

    :param job: Job to add the task to
    :param path: Path
    :param stem: Output file stem
    :param deps: Dependencies
    :param cpus: Number of CPUs
    :param mem: Amount of memory (GB)
    :param workers: How many workers to assign to this task
    :return: HyperQueue task
    """
    if workers > 1:
        log_paths = [
            log_path.with_stem(f"{log_path.stem}{i:02d}") for i in range(workers)
        ]
        deps = [
            _automech_hyperqueue_task(
                job=job,
                path=path,
                log_path=p,
                deps=deps,
                cpus=cpus,
                mem=mem,
                ignore_error=True,
            )
            for p in log_paths
        ]

    return _automech_hyperqueue_task(
        job=job, path=path, log_path=log_path, deps=deps, cpus=cpus, mem=mem
    )


def _automech_hyperqueue_task(
    job: Job,
    path: Path,
    log_path: Path,
    deps: Sequence[HQTask],
    cpus: int,
    mem: int,
    lock: bool = True,
    ignore_error: bool = False,
) -> HQTask:
    """Create a HyperQueue task to run automech."""
    run_ = run_automech

    # Create lock file if requested
    lock_path = log_path.with_suffix(Extension.running)
    run_ = lock_wrapper(run_, lock_file=lock_path) if lock else run_

    # Ignore errors if requested
    run_ = ignore_error_wrapper(run_) if ignore_error else run_

    return job.function(
        fn=run_,
        cwd=path,
        stdout=str(log_path),
        stderr=str(log_path),
        deps=deps,
        resources=ResourceRequest(cpus=cpus, resources={"mem": memory_mib(mem)}),
    )


# Ignore error
def ignore_error_wrapper(func: Callable[..., None]) -> Callable[..., None]:
    """Generate function wrapper that creates lock file during function execution.

    :param func: Function
    :param lock_file: Lock file to create while running
    :return: Wrapped function
    """

    @functools.wraps(func)
    def wrapper(*args, **kwargs) -> None:
        try:
            func(*args, **kwargs)
        except Exception as err:
            print(err)

    return wrapper


# Lock file (indicates that AutoMech is running)
def lock_wrapper(
    func: Callable[..., None], lock_file: str | Path
) -> Callable[..., None]:
    """Generate function wrapper that creates lock file during function execution.

    :param func: Function
    :param lock_file: Lock file to create while running
    :return: Wrapped function
    """
    lock_file = Path(lock_file)

    @functools.wraps(func)
    def wrapper(*args, **kwargs) -> None:
        with lock_file_context(lock_file):
            func(*args, **kwargs)

    return wrapper


@contextlib.contextmanager
def lock_file_context(lock_file: str | Path):
    lock_file = Path(lock_file)

    try:
        lock_file.touch()
        print(f"Created lock file {lock_file}")
        yield
    finally:
        lock_file.unlink()
        print(f"Removed lock file {lock_file}")


def dependency_graph(task_groups: Sequence[Sequence[Task]]) -> nx.DiGraph:
    """Create a subtask dependency graph from task groups."""
    # Store the task count and the subtasks keys for each group, for determining
    # group-level dependencies
    group_dct = {
        group_idx: (len(tasks) - 1, [s.key for s in tasks[0].subtasks])
        for group_idx, tasks in enumerate(task_groups)
    }
    group_idx0_dct = {
        group_idx: next(
            (i for i in reversed(range(group_idx)) if all(group_dct[i])), None
        )
        for group_idx, _ in enumerate(task_groups)
    }

    # Build the dependency graph
    dep_graph = nx.DiGraph()
    for group_idx, tasks in enumerate(task_groups):
        # Add dependencies within the group
        for task_idx, task in enumerate(tasks):
            dep_task_idx = task_idx - 1 if task_idx > 0 else None
            for subtask in task.subtasks:
                dep_graph.add_node((group_idx, task_idx, subtask.key))
                if dep_task_idx is not None:
                    dep_graph.add_edge(
                        (group_idx, dep_task_idx, subtask.key),
                        (group_idx, task_idx, subtask.key),
                    )

        # Add dependencies between groups
        group_idx0 = group_idx0_dct.get(group_idx)
        if group_idx0 is not None:
            task_idx0, subtask_keys0 = group_dct[group_idx0]
            _, subtask_keys = group_dct[group_idx]
            task_idx = 0
            for key0, key in itertools.product(subtask_keys0, subtask_keys):
                dep_graph.add_edge(
                    (group_idx0, task_idx0, key0), (group_idx, task_idx, key)
                )

    assert nx.is_weakly_connected(dep_graph), (
        "Dependency graph must not be disconnected:\n"
        f"group_idx0_dct = {group_idx0_dct}\ngroup_dct={group_dct}"
    )

    return dep_graph


def memory_mib(mem: int) -> int:
    """Convert memory in GB to MiB.

    :param mem: Memory (GB)
    :return: Memory (MiB)
    """
    return math.ceil(pint.Quantity(mem, "GB").m_as("MiB"))


def start_hyperqueue_server() -> None:
    """Re-start HyperQueue server."""
    print("Starting HyperQueue server...")
    subprocess.Popen(
        ["hq", "server", "start"], stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
    )
    # Wait up to 1 second for the file to appear
    for _ in range(10):
        time.sleep(0.1)
        if os.path.exists(HQ_PATH):
            break
    assert os.path.exists(HQ_PATH), f"Could not start server at {HQ_PATH}"


def add_hyperqueue_allocation(mem: int, cpus: int, flags: str) -> None:
    """Create a HyperQueue allocation.

    :param mem: Memory (GB)
    :param nprocs: Number of processers
    :param flags: Additional flags for sbatch/qsub
    """
    print(f"Adding HyperQueue allocation with mem={mem}GB, cpus={cpus}, flags={flags}")

    if shutil.which("sbatch"):
        print("Detected SLURM on system. HyperQueue allocation command:")
        args = [
            "hq",
            "alloc",
            "add",
            "slurm",
            "--time-limit",
            "1h",
            f"--cpus={cpus}",
            f"--resource=mem=sum({memory_mib(mem)})",
            "--",
            f"--mem={mem}G",
            "--ntasks=1",
            *flags.split(),
        ]
        print(" ".join(args))
        subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    elif shutil.which("qsub"):
        print("Detected PBS on system. HyperQueue allocation command:")
        args = [
            "hq",
            "alloc",
            "add",
            "pbs",
            "--time-limit",
            "1h",
            f"--cpus={cpus}",
            f"--resource=mem=sum({memory_mib(mem)})",
            "--",
            *flags.split(),
        ]
        print(" ".join(args))
        subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        msg = (
            "No SLURM or PBS detected. Please manually configure HyperQueue allocation."
        )
        raise ValueError(msg)
