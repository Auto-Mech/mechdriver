""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

import itertools
from collections.abc import Sequence
from pathlib import Path

import more_itertools as mit
import pandas
import yaml

from ._check_log import STATUS_WIDTH, Status, colored_status_string, parse_log_status
from ._subtasks_setup import (
    INFO_FILE,
    SUBTASK_DIR,
    InfoKey,
    TableKey,
    Task,
    read_task_list,
)


def main(
    path: str | Path = SUBTASK_DIR,
) -> None:
    """Runs subtasks in parallel on Ad Hoc cluster

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param nodes: A comma-separated list of nodes to run on
    :param activation_hook: Shell commands for activating the AutoMech environment on the remote
    """
    path = Path(path)
    assert (
        path.exists()
    ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    info_path = path / INFO_FILE
    info_dct = yaml.safe_load(info_path.read_text())

    group_ids = info_dct[InfoKey.group_ids]
    work_path = info_dct[InfoKey.work_path]
    run_path = Path(info_dct[InfoKey.run_path])
    save_path = Path(info_dct[InfoKey.save_path])

    run_path.mkdir(exist_ok=True)
    save_path.mkdir(exist_ok=True)

    non_okay_log_records = []
    for group_id in group_ids:
        df = pandas.read_csv(path / f"{group_id}.csv")
        tasks = read_task_list(path / f"{group_id}.yaml")
        twidth = task_column_width(tasks)
        skeys = subtask_keys(tasks)

        header = task_group_header(skeys, twidth)
        print(header)
        for task_key, row in df.iterrows():
            task: Task = tasks[task_key]
            assert row[TableKey.task] == task.name, f"{row} does not match {task.name}"

            subtask_paths = list(map(row.get, skeys))
            subtask_abs_paths = [Path(work_path) / p for p in subtask_paths]
            subtask_stats = list(
                map(colored_status_string, map(parse_subtask_status, subtask_abs_paths))
            )
            line = f"{task.name:>{twidth}} " + " ".join(subtask_stats)
            print(line)

            non_okay_log_records.extend(
                (task.name, skey, p, s)
                for skey, spath in zip(skeys, subtask_paths, strict=True)
                for p, s in log_paths_with_statuses(
                    spath, exclude_stats=[Status.OK]
                ).items()
            )
        print()

    if non_okay_log_records:
        print(f"Non-OK log files in {work_path}:")
        twidth = max(len(r[0]) for r in non_okay_log_records)
        swidth = max(len(r[1]) for r in non_okay_log_records)
        pwidth = max(len(r[2]) for r in non_okay_log_records)
        for task_name, skey, log_path, stat in non_okay_log_records:
            stat = colored_status_string(stat)
            print(f"{task_name:<{twidth}} {skey:<{swidth}} {log_path:<{pwidth}} {stat}")
    print()


def log_paths_with_statuses(
    path: str | Path, exclude_stats: Sequence[Status] = ()
) -> dict[str, Status]:
    """Get a dictionary of log file paths and statuses at a given path

    :param path: _description_
    :return: _description_
    """
    log_paths = list(map(str, Path(path).glob("out*.log")))
    log_stats = list(map(parse_log_status, log_paths))
    log_dct = dict(zip(log_paths, log_stats, strict=True))
    return {k: v for k, v in log_dct.items() if v not in exclude_stats}


def parse_subtask_status(path: str | Path, small_thresh: float = 0.2) -> Status:
    """Parse the run status from a subtask directory

    :param path: The directory path
    :return: The status
    """
    log_dct = log_paths_with_statuses(path)
    if not log_dct:
        return Status.TBD

    log_stats = list(log_dct.values())
    log_stat_set = set(log_stats)

    # All log files have the same status -> <common status>
    if len(log_stat_set) == 1:
        return next(iter(log_stat_set))

    # Some log files are still runnning -> RUNNING
    if Status.RUNNING in log_stat_set:
        return Status.RUNNING

    # Some log files have errors -> ERROR | OKAY_1E | OKAY_2E
    error_count = log_stats.count(Status.ERROR)
    error_frac = error_count / len(log_stats)
    if error_count == 1 and error_frac < small_thresh:
        return Status.OK_1E
    if error_count == 2 and error_frac < small_thresh:
        return Status.OK_2E
    if Status.ERROR in log_stat_set:
        return Status.ERROR

    # Some log fils have warnings -> WARNING
    assert log_stat_set == {Status.OK, Status.WARNING}
    return Status.WARNING


def task_column_width(tasks: list[Task]) -> int:
    """Get the appropriate column width for a list of tasks

    :param tasks: The list of tasks
    :return: The column width
    """
    return max(map(len, (task.name for task in tasks)))


def subtask_keys(tasks: list[Task]) -> list[str]:
    """Get the list of subtask keys

    If tasks have different sets of subtask keys, this returns the union of all of them

    :param tasks: The list of tasks
    :return: The subtask keys
    """
    return list(mit.unique_everseen(itertools.chain(*(t.subtask_keys for t in tasks))))


def task_group_header(skeys: list[str], twidth: int) -> str:
    """Print the header for a task group

    :param tasks: The list of tasks
    :param twidth: The task column width
    :return: The header
    """
    header = f"{TableKey.task:>{twidth}} "
    header += " ".join(f"{k:^{STATUS_WIDTH}}" for k in skeys)
    return header
