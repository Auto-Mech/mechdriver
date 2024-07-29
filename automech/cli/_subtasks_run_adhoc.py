""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

import subprocess
from pathlib import Path

import pandas
import yaml

from ._subtasks_setup import INFO_FILE, SUBTASK_DIR, InfoKey, SpecKey, TableKey

SCRIPT_DIR = Path(__file__).parent / "scripts"
RUN_SCRIPT = str(SCRIPT_DIR / "run_adhoc.sh")


def main(
    path: str | Path = SUBTASK_DIR,
    nodes: str | None = None,
    activation_hook: str | None = None,
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

    for group_id in group_ids[:1]:
        df = pandas.read_csv(path / f"{group_id}.csv")
        spec_lst = yaml.safe_load((path / f"{group_id}.yaml").read_text())
        for task_key, row in df.iterrows():
            task_name = row[TableKey.task]
            spec_dct: dict = spec_lst[task_key]
            assert spec_dct[SpecKey.task] == task_name, f"{task_name} {spec_dct}"

            subtask_nprocs = str(spec_dct.get(SpecKey.nprocs))
            subtask_mem = str(spec_dct.get(SpecKey.mem))
            subtask_paths = ",".join(v for k, v in row.items() if str(k).isdigit())

            run_args = [
                RUN_SCRIPT,
                work_path,
                subtask_mem,
                subtask_nprocs,
                subtask_paths,
                nodes,
                "" if activation_hook is None else activation_hook,
            ]
            subprocess.run(run_args)
