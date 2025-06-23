"""Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster"""

import itertools
import os
import subprocess
from collections.abc import Sequence
from pathlib import Path

import networkx as nx
import pandas as pd
import yaml
from hyperqueue import Client, Job
from hyperqueue.ffi.protocol import ResourceRequest

from ..base import Status
from ._0setup import INFO_FILE, SUBTASK_DIR, SubtasksInfo, Task
from ._1status import log_paths_with_check_results, parse_subtask_status

SCRIPT_DIR = Path(__file__).parent / "scripts"
RUN_SCRIPT = str(SCRIPT_DIR / "run_adhoc.sh")
HOME = Path(os.environ.get("HOME"))


def run(
    path: str = ".",
    dir_name: str = SUBTASK_DIR,
    hyperqueue_path: str | None = None,
    statuses: Sequence[Status] = (Status.TBD,),
) -> None:
    """Run subtasks in parallel using HyperQueue.

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param hyperqueue_path: The path to the HyperQueue server directory
    :param statuses: A comma-separated list of status to run or re-run
    :param tar: Tar the subtask data and save filesystem after running?
    """
    path = Path(path).resolve()
    dir_path = path / dir_name
    info_file = dir_path / INFO_FILE
    info = SubtasksInfo.model_validate(yaml.safe_load(info_file.read_text()))

    # Set up the HyperQueue client
    hyperqueue_path = hyperqueue_path or HOME / ".hq-server" / "hq-current"
    client = Client(hyperqueue_path)

    # Set up the HyperQueue job workflow
    job = Job()

    # For now, just do this for the first task group
    job_dct = {}
    dep_graph = dependency_graph(info.task_groups)
    for group_idx, task_group in enumerate(info.task_groups):
        for task_idx, task in enumerate(task_group):
            for subtask in task.subtasks:
                # Determine dependencies from dependency graph
                job_key = (group_idx, task_idx, subtask.key)
                dep_job_keys = dep_graph.predecessors(job_key)
                deps = list(map(job_dct.get, dep_job_keys))

                # Create the job
                subtask_path = dir_path / subtask.path
                stem = "out"
                subtask_job = job.program(
                    ["automech", "run", "-p", str(subtask_path), "-r", stem],
                    cwd=subtask_path,
                    stdout=subtask_path / f"{stem}.log",
                    stderr=subtask_path / f"{stem}.log",
                    deps=deps,
                    resources=ResourceRequest(resources={"mem": task.mem * 1000}),
                )

                # Add the job to the job dictionary
                job_dct[job_key] = subtask_job

    submitted = client.submit(job)
    client.wait_for_jobs([submitted])


def dependency_graph(task_groups: Sequence[Sequence[Task]]) -> nx.DiGraph:
    """Create a subtask dependency graph from task groups."""
    # Turn task groups into dataframes, and add a "keys" column
    task_dfs = [pd.DataFrame([t.model_dump() for t in ts]) for ts in task_groups]
    for task_df in task_dfs:
        task_df["keys"] = task_df["subtasks"].apply(
            lambda subtasks: [s["key"] for s in subtasks]
        )

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
        if group_idx > 0:
            group_idx0 = group_idx - 1
            task_df0 = task_dfs[group_idx0]
            task_idx0 = task_df0.index[-1]
            subtask_keys0 = task_df0.iloc[task_idx0]["keys"]

            task_df = task_dfs[group_idx]
            task_idx = 0
            subtask_keys = task_df.iloc[task_idx]["keys"]
            for key0, key in itertools.product(subtask_keys0, subtask_keys):
                dep_graph.add_edge(
                    (group_idx0, task_idx0, key0), (group_idx, task_idx, key)
                )

    return dep_graph


def run_multiple(
    paths: Sequence[str | Path] = (".",),
    nodes: Sequence[str] | None = None,
    dir_name: str = SUBTASK_DIR,
    activation_hook: str | None = None,
    statuses: Sequence[Status] = (Status.TBD,),
) -> None:
    """Runs multiple sets of subtasks in parallel on Ad Hoc cluster

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param paths: The paths where the AutoMech subtasks were set up
    :param nodes: A list of nodes to run on
    :param activation_hook: Shell commands for activating the AutoMech environment on the remote
    :param statuses: A comma-separated list of status to run or re-run
    :param tar: Tar the subtask data and save filesystem after running?
    """
    # Determine paths
    paths = [Path(p).resolve() for p in paths]
    dir_paths = [p / dir_name for p in paths]
    info_files = [d / INFO_FILE for d in dir_paths]
    for dir_path in dir_paths:
        assert (
            dir_path.exists()
        ), f"Path not found: {dir_path}.\nDid you run `automech subtasks setup` first?"

    # Read in subtask information
    infos = [SubtasksInfo(**yaml.safe_load(f.read_text())) for f in info_files]

    # Make sure the run and save directories exist
    for path, info in zip(paths, infos, strict=True):
        (path / info.run_path).mkdir(exist_ok=True)
        (path / info.save_path).mkdir(exist_ok=True)

    # Zip tasks together in sequence
    tasks_lst = list(
        itertools.zip_longest(*(itertools.chain(*i.task_groups) for i in infos))
    )

    for tasks in tasks_lst:
        mem = 0
        nprocs = 0
        work_paths = []
        subtask_paths = []
        subtask_logs = []
        for path, task in (
            (p, t) for p, t in zip(paths, tasks, strict=True) if t is not None
        ):
            mem = max(mem, task.mem)
            nprocs = max(nprocs, task.nprocs)
            for subtask in task.subtasks:
                subtask_path = dir_name / subtask.path
                subtask_path_abs = path / subtask_path
                status = parse_subtask_status(
                    log_paths_with_check_results(subtask_path_abs)
                )
                if status in statuses:
                    work_paths.extend([path] * subtask.nworkers)
                    subtask_paths.extend([subtask_path] * subtask.nworkers)
                    subtask_logs.extend(
                        [f"out{i}.log" for i in range(subtask.nworkers)]
                    )

        if subtask_paths:
            run_args = [
                RUN_SCRIPT,
                ",".join(map(str, work_paths)),
                f"{mem}",
                f"{nprocs}",
                ",".join(map(str, subtask_paths)),
                ",".join(subtask_logs),
                ",".join(nodes),
                "" if activation_hook is None else activation_hook,
            ]
            subprocess.run(run_args)
