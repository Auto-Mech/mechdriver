""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

import itertools
from collections.abc import Sequence
from pathlib import Path

import more_itertools as mit
import yaml

from ..base import Status, check_log, colored_status_string
from ._0setup import INFO_FILE, SUBTASK_DIR, SubtasksInfo, Task


def status(
    path: str | Path = ".",
    dir_name: str = SUBTASK_DIR,
    check_file: str = "check.log",
    wrap: int = 18,
) -> None:
    """Check the status of running subtasks

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param check_file: Log file for writing paths to be checked
    :Param wrap: Wrap to include this many subtask columns per row
    """
    path = Path(path).resolve()
    assert (
        path.exists()
    ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    info_file = path / dir_name / INFO_FILE
    info = SubtasksInfo(**yaml.safe_load(info_file.read_text()))

    check_records = []
    for tasks in info.task_groups:
        skeys = subtask_keys(tasks)
        kwidth = key_column_width(tasks)
        vwidth = value_column_width(tasks)

        print_long_row_guide(
            kwidth=kwidth, vwidth=vwidth, nvals=len(skeys), wrap=wrap, char="#"
        )
        print_task_row("task", skeys, kwidth=kwidth, vwidth=vwidth, wrap=wrap)
        for task in tasks:
            vals = []
            for subtask in task.subtasks:
                log_dct = log_paths_with_check_results(path / dir_name / subtask.path)
                stat = parse_subtask_status(log_dct=log_dct)
                vals.append(colored_status_string(stat, width=vwidth))
                check_records.extend(
                    (task.name, subtask.key, p, s, L)
                    for p, (s, L) in log_dct.items()
                    if s != Status.OK
                )
            print_task_row(task.name, vals, kwidth=kwidth, vwidth=vwidth, wrap=wrap)
        print()

    check_lines = []
    if check_records:
        check_lines.append(f"Non-OK log files in {path / dir_name}:")
        kwidth = max(len(r[0]) for r in check_records)
        vwidth = max(len(r[1]) for r in check_records)
        pwidth = max(len(r[2]) for r in check_records)
        for task_name, skey, log_path, stat, line in check_records:
            stat = colored_status_string(stat)
            check_lines.append(
                f"~{task_name:<{kwidth}} {skey:<{vwidth}} {log_path:<{pwidth}} {stat}"
            )

            if line is not None:
                check_lines.append(line)

    check_file: Path = Path(check_file)
    check_file_contents = "\n".join(check_lines)
    check_file.write_text(f"{check_file_contents}\n" if check_file_contents else "")


def log_paths_with_check_results(
    path: str | Path,
) -> dict[str, tuple[Status, str | None]]:
    """Get a dictionary of log file paths and statuses at a given path

    :param path: The directory path
    :return: A dictionary mapping log paths onto log check results
    """
    log_paths = list(map(str, Path(path).glob("out*.log")))
    if not log_paths:
        return {}

    log_checks = list(map(check_log, log_paths))
    return dict(zip(log_paths, log_checks, strict=True))


def parse_subtask_status(
    log_dct: dict[str, tuple[Status, str | None]], small_thresh: float = 0.2
) -> Status:
    """Parse the run status from a subtask directory

    :return: The status
    """
    if not log_dct:
        return Status.TBD

    log_stats, *_ = zip(*log_dct.values())
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


def subtask_keys(tasks: list[Task]) -> list[str]:
    """Get the list of subtask keys

    If tasks have different sets of subtask keys, this returns the union of all of them

    :param tasks: The list of tasks
    :return: The subtask keys
    """
    return list(
        mit.unique_everseen(
            itertools.chain(*([s.key for s in t.subtasks] for t in tasks))
        )
    )


def key_column_width(tasks: list[Task]) -> int:
    """Get the appropriate column width for a list of tasks

    :param tasks: The list of tasks
    :return: The column width
    """
    return max(map(len, (task.name for task in tasks)))


def value_column_width(tasks: list[Task]) -> int:
    """Get the appropriate column width for a list of tasks

    :param tasks: The list of tasks
    :return: The column width
    """
    skeys = subtask_keys(tasks)
    return max(*map(len, skeys), *(len(s.value) for s in Status))


def print_task_row(
    key: str, vals: Sequence[str], kwidth: int, vwidth: int, wrap: int
) -> None:
    """Print a single row in the task group table

    :param key: The row label
    :param vals: The row values
    :param kwidth: The label column width
    :param vwidth: The value column width
    :param wrap: Wrap the row values after this many columns
    """
    for chunk_vals in mit.chunked(vals, wrap):
        row = f"{key:>{kwidth}} "
        row += " ".join(f"{v:^{vwidth}}" for v in chunk_vals)
        print(row)
        key = ""  # drop the label after the first chunk

    # If wrapping, add an extra dividing line as a guide
    print_long_row_guide(kwidth=kwidth, vwidth=vwidth, nvals=len(vals), wrap=wrap)


def print_long_row_guide(
    kwidth: int, vwidth: int, nvals: int, wrap: int, char: str = "-"
) -> None:
    """Print a horizontal guide to guide the eye, if the row is long

    :param kwidth: The label column width
    :param vwidth: The value column width
    :param wrap: Wrap the row values after this many columns
    :param char: The character to use for the separator, defaults to "-"
    """
    if nvals > wrap:
        total_width = kwidth + 1 + (vwidth + 1) * wrap
        print(char * total_width)
