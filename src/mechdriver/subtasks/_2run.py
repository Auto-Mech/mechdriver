"""Standalone script to run AutoMech subtasks in parallel using HyperQueue."""

import contextlib
import functools
import itertools
from collections.abc import Callable, Sequence
from pathlib import Path

import networkx as nx
import yaml

from ..base import Extension, Status
from ..base import run as run_automech
from . import hq
from ._0setup import INFO_FILE, SUBTASK_DIR, SubtasksInfo, Task
from ._1status import log_paths_with_check_results, parse_subtask_status

FunctionKey = tuple[int, int, str]


def run_multiple(
    paths: Sequence[str | Path] = (".",),
    dir_name: str = SUBTASK_DIR,
    statuses: Sequence[Status] = (Status.TBD,),
    manager: str | None = None,
    manager_flags: str | None = None,
    env_prologue: str | None = None,
) -> None:
    """Run multiple sets of subtasks in parallel using HyperQueue.

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param paths: The paths where the AutoMech subtasks were set up
    :param dir_name: The subtask directory name
    :param hyperqueue_path: The path to the HyperQueue server directory
    :param statuses: A comma-separated list of status to run or re-run
    :param manager: Specify workload manager (PBS or Slurm) instead of autodetecting
    :param manager_flags: Automatically configure HyperQueue allocation with
        these PBS or Slurm submission flags
    :param env_prologue: Command(s) to activate Python environment
    """
    # If flags were passed in, attempt to auto-configure
    if manager_flags is not None:
        # Start HyperQueue server
        hq.start_server()

        # Determine max memory and CPU requirements across all paths
        grouped_tasks = [subtasks_info_tasks(p, dir_name=dir_name) for p in paths]
        flat_tasks = list(itertools.chain.from_iterable(grouped_tasks))
        mem = max(t.mem for t in flat_tasks)
        cpus = max(t.nprocs for t in flat_tasks)

        # Start auto-allocation queue
        hq.create_allocation_queue(
            mem=mem, cpus=cpus, flags=manager_flags, manager=manager
        )

    # Create HyperQueue client
    client = hq.client(env_prologue=env_prologue)

    # Set up the HyperQueue job workflow
    job = hq.job()

    for path in paths:
        job = setup_job(path=path, dir_name=dir_name, statuses=statuses, job=job)

    submitted_job = client.submit(job)
    client.wait_for_jobs([submitted_job])


def setup_job(
    path: str | Path = ".",
    dir_name: str = SUBTASK_DIR,
    statuses: Sequence[Status] = (Status.TBD,),
    job: hq.Job | None = None,
) -> hq.Job:
    """Run subtasks in parallel using HyperQueue.

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param dir_name: The subtask directory name
    :param statuses: A comma-separated list of status to run or re-run
    :param job: Append to an existing job
    """
    sub_path = subtasks_path(path=path, dir_name=dir_name)
    sub_info = subtasks_info(path=path, dir_name=dir_name)

    # Make sure the run and save directories exist
    sub_info.run_path.mkdir(exist_ok=True)
    sub_info.save_path.mkdir(exist_ok=True)

    # Set up the HyperQueue job workflow
    job = job or hq.job()

    # For now, just do this for the first task group
    func_dct: dict[FunctionKey, hq.Function] = {}
    dep_graph = dependency_graph(sub_info.task_groups)
    for group_idx, task_group in enumerate(sub_info.task_groups):
        for task_idx, task in enumerate(task_group):
            for subtask in task.subtasks:
                # Determine dependencies from dependency graph
                func_key = (group_idx, task_idx, subtask.key)
                dep_func_keys = dep_graph.predecessors(func_key)
                dep_funcs = [func_dct[k] for k in dep_func_keys]

                func = assign_function(
                    job=job,
                    path=sub_path / subtask.path,
                    log_path=sub_path / subtask.path / "out.log",
                    deps=dep_funcs,
                    cpus=task.nprocs,
                    mem=task.mem,
                    workers=subtask.nworkers,
                )

                # Add the job to the job dictionary
                func_dct[func_key] = func

    return job


def assign_function(
    job: hq.Job,
    path: Path,
    log_path: Path,
    deps: Sequence[hq.Function],
    cpus: int,
    mem: int,
    workers: int = 1,
) -> hq.Function:
    r"""Assign function(s) to HyperQueue job.

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
            assign_atomic_function(
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

    return assign_atomic_function(
        job=job, path=path, log_path=log_path, deps=deps, cpus=cpus, mem=mem
    )


def assign_atomic_function(
    job: hq.Job,
    path: Path,
    log_path: Path,
    deps: Sequence[hq.Function],
    cpus: int,
    mem: int,
    lock: bool = True,
    ignore_error: bool = False,
) -> hq.Function:
    """Create a HyperQueue task to run automech."""
    run_ = run_automech

    # Create lock file if requested
    lock_path = log_path.with_suffix(Extension.running)
    run_ = lock_wrapper(run_, lock_file=lock_path) if lock else run_

    # Ignore errors if requested
    run_ = ignore_error_wrapper(run_) if ignore_error else run_

    resources = hq.resource_request(cpus=cpus, mem=mem)
    stdout = stderr = str(log_path)
    return job.function(
        fn=run_, cwd=path, stdout=stdout, stderr=stderr, deps=deps, resources=resources
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


# Helpers
def subtasks_path(path: str | Path = ".", dir_name: str = SUBTASK_DIR) -> Path:
    """Determine absolute path to subtasks directory."""
    return Path(path).resolve() / dir_name


def subtasks_info(path: str | Path = ".", dir_name: str = SUBTASK_DIR) -> SubtasksInfo:
    """Read subtasks info from subtasks directory."""
    sub_info_path = subtasks_path(path=path, dir_name=dir_name) / INFO_FILE
    return SubtasksInfo.model_validate(yaml.safe_load(sub_info_path.read_text()))


def subtasks_info_tasks(
    path: str | Path = ".", dir_name: str = SUBTASK_DIR
) -> list[Task]:
    """Read subtasks info from subtasks directory."""
    sub_info = subtasks_info(path=path, dir_name=dir_name)
    return list(itertools.chain.from_iterable(sub_info.task_groups))
