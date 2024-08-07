""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

import subprocess
from pathlib import Path

import pandas
import yaml

from ._check_log import Status
from ._subtasks_setup import (
    INFO_FILE,
    SUBTASK_DIR,
    InfoKey,
    TableKey,
    Task,
    read_task_list,
)
from ._subtasks_status import parse_subtask_status

SCRIPT_DIR = Path(__file__).parent / "scripts"
RUN_SCRIPT = str(SCRIPT_DIR / "run_adhoc.sh")


def main(
    path: str | Path = SUBTASK_DIR,
    nodes: str | None = None,
    activation_hook: str | None = None,
    statuses: str = f"{Status.TBD.value}",
) -> None:
    """Runs subtasks in parallel on Ad Hoc cluster

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    :param nodes: A comma-separated list of nodes to run on
    :param activation_hook: Shell commands for activating the AutoMech environment on the remote
    :param statuses: A comma-separated list of status to run or re-run
    """
    path = Path(path)
    assert (
        path.exists()
    ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    statuses = list(map(Status, statuses.split(",")))

    info_path = path / INFO_FILE
    info_dct = yaml.safe_load(info_path.read_text())

    group_ids = info_dct[InfoKey.group_ids]
    work_path = info_dct[InfoKey.work_path]
    run_path = Path(info_dct[InfoKey.run_path])
    save_path = Path(info_dct[InfoKey.save_path])

    run_path.mkdir(exist_ok=True)
    save_path.mkdir(exist_ok=True)

    for group_id in group_ids:
        df = pandas.read_csv(path / f"{group_id}.csv")
        tasks = read_task_list(path / f"{group_id}.yaml")
        for task_key, row in df.iterrows():
            task: Task = tasks[task_key]
            assert row[TableKey.task] == task.name, f"{row} does not match {task.name}"

            subtask_paths = []
            subtask_logs = []

            for key, nworkers in zip(
                task.subtask_keys, task.subtask_nworkers, strict=True
            ):
                assert key in row, f"Key {key} not present in row:\n{row}"
                subtask_path = row.get(key)
                status = parse_subtask_status(subtask_path)
                if status in statuses:
                    subtask_paths.extend([subtask_path] * nworkers)
                    subtask_logs.extend([f"out{i}.log" for i in range(nworkers)])

            if subtask_paths:
                run_args = [
                    RUN_SCRIPT,
                    work_path,
                    f"{task.mem}",
                    f"{task.nprocs}",
                    ",".join(subtask_paths),
                    ",".join(subtask_logs),
                    nodes,
                    "" if activation_hook is None else activation_hook,
                ]
                subprocess.run(run_args)
