""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

import itertools
import subprocess
from collections.abc import Sequence
from pathlib import Path

import yaml

from ..base import Status
from ._0setup import INFO_FILE, SUBTASK_DIR, SubtasksInfo
from ._1status import log_paths_with_check_results, parse_subtask_status

SCRIPT_DIR = Path(__file__).parent / "scripts"
RUN_SCRIPT = str(SCRIPT_DIR / "run_adhoc.sh")


def run(
    path: str = ".",
    nodes: Sequence[str] | None = None,
    dir_name: str = SUBTASK_DIR,
    activation_hook: str | None = None,
    statuses: Sequence[Status] = (Status.TBD,),
) -> None:
    """Runs subtasks in parallel on Ad Hoc cluster

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param nodes: A list of nodes to run on
    :param activation_hook: Shell commands for activating the AutoMech environment on the remote
    :param statuses: A comma-separated list of status to run or re-run
    :param tar: Tar the subtask data and save filesystem after running?
    """
    return run_multiple(
        paths=[path],
        nodes=nodes,
        dir_name=dir_name,
        activation_hook=activation_hook,
        statuses=statuses,
    )


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
                status = parse_subtask_status(
                    log_paths_with_check_results(subtask_path)
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
