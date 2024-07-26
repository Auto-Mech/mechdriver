""" Standalone script to run AutoMech subtasks in parallel on an Ad Hoc SSH Cluster
"""

from pathlib import Path

import pandas
import yaml

from ._subtasks_setup import SUBTASK_DIR, INFO_FILE, InfoKey


def main(path: str | Path=SUBTASK_DIR) -> None:
    """Runs subtasks in parallel on Ad Hoc cluster

    Assumes the subtasks were set up at this path using `automech subtasks setup`

    :param path: The path where the AutoMech subtasks were set up
    """
    path = Path(path)
    assert (
        path.exists()
    ), f"Path not found: {path}.\nDid you run `automech subtasks setup` first?"

    info_path = path / INFO_FILE
    info_dct = yaml.safe_load(info_path.read_text())

    group_ids = info_dct[InfoKey.group_ids]

    print(info_dct)
    print(group_ids)

    for group_id in group_ids[:1]:
        group_csv = path / f"{group_id}.csv"
        group_df = pandas.read_csv(group_csv)
        print(group_df)
