"""Standalone script to break an AutoMech input into subtasks for parallel execution"""

from collections.abc import Sequence
from pathlib import Path

import pydantic
import yaml

from . import util

DEFAULT_TASK_GROUPS = ("els", "thermo", "ktp", "proc")
GROUP_ID = {"els-spc": 0, "els-pes": 1, "thermo": 2, "ktp": 3, "proc-spc": 4, "proc-pes": 5}
GROUP_TASK_AND_KEY_TYPE = {
    "els-spc": ("els", "spc"),
    "els-pes": ("els", "pes"),
    "thermo": ("thermo", "spc"),
    "ktp": ("ktp", "pes"),
    "proc-spc": ("proc", "spc"),
    "proc-pes": ("proc", "pes"),
}
COMBINED_TASK_GROUPS = ()#"thermo", "ktp")

SUBTASK_DIR = "subtasks"
INFO_FILE = "info.yaml"


class Subtask(pydantic.BaseModel):
    key: str
    nworkers: int
    path: Path

    @pydantic.field_serializer("path")
    def serialize_path(self, path: Path, _):
        return str(path)

    def __rtruediv__(self, path: str | Path):
        return self.model_copy(update={"path": path / self.path})


class Task(pydantic.BaseModel):
    name: str
    line: str
    mem: int
    nprocs: int
    subtasks: list[Subtask]

    def __rtruediv__(self, path: str | Path):
        return self.model_copy(update={"subtasks": [path / s for s in self.subtasks]})


class SubtasksInfo(pydantic.BaseModel):
    run_path: Path
    save_path: Path
    task_groups: list[list[Task]]

    @pydantic.field_serializer("run_path", "save_path")
    def serialize_path(self, path: Path):
        return str(path)

    def __rtruediv__(self, path: str | Path):
        return self.model_copy(
            update={"task_groups": [[path / t for t in ts] for ts in self.task_groups]}
        )


def setup_multiple(
    paths: Sequence[str | Path],
    dir_name: str = SUBTASK_DIR,
    save_path: str | Path | None = None,
    run_path: str | Path | None = None,
    task_group_keys: Sequence[str] = DEFAULT_TASK_GROUPS,
):
    """Creates run directories for each task/species/TS and returns the paths in tables

    Task types: 'els', 'thermo', or 'ktp'
    Subtask types: 'spc', 'pes', or `None` (=all species and/or reactions)

    :param paths: The paths to the AutoMech inputs to split into subtasks
    :param dir_name: The name of the subtask directory to write to
        (will be created at `path`)
    :param save_path: The path to the save filesystem
        (if `None`, the value in run.dat is used)
    :param run_path: The path to the run filesystem
        (if `None`, the value in run.dat is used)
    :param task_groups: The task groups to set up
    :return: DataFrames of run paths, whose columns (species/TS index) are independent
        and can be run in parallel, but whose rows (tasks) are potentially sequential
    """
    for path in paths:
        setup(
            path=path,
            dir_name=dir_name,
            save_path=save_path,
            run_path=run_path,
            task_group_keys=task_group_keys,
        )


def setup(
    path: str | Path,
    dir_name: str = SUBTASK_DIR,
    save_path: str | Path | None = None,
    run_path: str | Path | None = None,
    absolute: bool = True,
    task_group_keys: Sequence[str] = DEFAULT_TASK_GROUPS,
):
    """Creates run directories for each task/species/TS and returns the paths in tables

    Task types: 'els', 'thermo', or 'ktp'
    Subtask types: 'spc', 'pes', or `None` (=all species and/or reactions)

    :param path: The path to the AutoMech input to split into subtasks
    :param dir_name: The name of the subtask directory to write to
        (will be created at `path`)
    :param save_path: The path to the save filesystem
        (if `None`, the value in run.dat is used)
    :param run_path: The path to the run filesystem
        (if `None`, the value in run.dat is used)
    :param absolute: Whether to resolve absolute paths
    :param task_groups: The task groups to set up
    :return: DataFrames of run paths, whose columns (species/TS index) are independent
        and can be run in parallel, but whose rows (tasks) are potentially sequential
    """
    # Pre-process the task groups
    task_group_keys = list(task_group_keys)
    if "els" in task_group_keys:
        idx = task_group_keys.index("els")
        task_group_keys[idx : idx + 1] = ("els-spc", "els-pes")
    if "proc" in task_group_keys:
        idx = task_group_keys.index("proc")
        task_group_keys[idx : idx + 1] = ("proc-spc", "proc-pes")
    assert all(g in GROUP_ID for g in task_group_keys), (
        f"{task_group_keys} not in {GROUP_ID}"
    )

    # Read input files from source path
    path = Path(path)
    print("----")
    print(f"Setting up subtasks at path {path.resolve()}")

    file_dct = util.read_input_files(path)
    run_dct = util.parse_run_dat(file_dct.get("run.dat"))

    # Set the run and save paths
    inp_dct = util.input_arguments_from_run_dict(
        run_dct, save_path=save_path, run_path=run_path, absolute=absolute
    )
    run_dct["input"] = "\n".join(f"{k} = {v}" for k, v in inp_dct.items())

    # Create the path for the subtask directories
    dir_path = path / dir_name
    dir_path.mkdir(exist_ok=True)

    # Set up the subtasks for each group
    task_groups = []
    for task_group_key in task_group_keys:
        group_id = GROUP_ID.get(task_group_key)
        task_type, key_type = GROUP_TASK_AND_KEY_TYPE.get(task_group_key)
        tasks = setup_subtask_group(
            run_dct,
            file_dct,
            task_type=task_type,
            key_type=key_type,
            group_id=group_id,
            dir_path=dir_path,
        )
        if tasks is not None:
            task_groups.append(tasks)

    # Write the subtask info to YAML
    info_path = dir_path / INFO_FILE
    print(f"Writing subtask information to {info_path}")
    print()
    info = SubtasksInfo(
        save_path=inp_dct.get("save_prefix"),
        run_path=inp_dct.get("run_prefix"),
        task_groups=task_groups,
    )
    info_path.write_text(yaml.safe_dump(info.model_dump()))


def setup_subtask_group(
    run_dct: dict[str, str],
    file_dct: dict[str, str],
    task_type: str,
    key_type: str | None = None,
    group_id: str | int | None = None,
    dir_path: str | Path = SUBTASK_DIR,
) -> list[Task] | None:
    """Set up a group of subtasks from a run dictionary, creating run directories and
    returning them in a table

    :param source_path: The path to the AutoMech input to split into subtasks
    :param task_type: The type of task: 'els', 'thermo', or 'ktp'
    :param key_type: The type of subtask key: 'spc', 'pes', or `None`
    :param group_id: The group ID, used to name files and folders
    :return: A DataFrame of run paths, whose columns (subtasks) are independent and can
        be run in parallel, but whose rows (tasks) are potentially sequential
    """
    # Form a prefix for the task/subtask type
    if group_id is None:
        type_keys = [task_type] + ([] if key_type is None else [key_type])
        return "_".join(type_keys)
    group_id = str(group_id)

    # Blocks that must be included in the run.dat
    block_keys = ["input"] + (["pes", "spc"] if key_type is None else [key_type])

    # Determine the task list
    tasks = determine_task_list(
        run_dct,
        file_dct,
        task_type,
        key_type=key_type,
        group_id=group_id,
    )

    # If the task list is empty, return `None`
    if not tasks:
        return None

    # Create directories for each subtask and save the paths in a DataFrame
    for task in tasks:
        print(f"Setting up subtask directories for {task.name}")
        task_run_dct = {k: v for k, v in run_dct.items() if k in block_keys}
        task_run_dct[task_type] = task.line
        for subtask in task.subtasks:
            # Generate the subtask path
            subtask_path = dir_path / subtask.path
            subtask_path.mkdir(parents=True, exist_ok=True)
            # Generate the input file dictionary
            subtask_run_dct = task_run_dct.copy()
            if subtask.key != util.ALL_KEY:
                subtask_run_dct[key_type] = subtask.key
            subtask_file_dct = {
                **file_dct,
                "run.dat": util.form_run_dat(subtask_run_dct),
            }
            # Write the input files and append the path to the current dataframe row
            util.write_input_files(subtask_path, subtask_file_dct)

    return tasks


def determine_task_list(
    run_dct: dict[str, str],
    file_dct: dict[str, str],
    task_type: str,
    key_type: str | None = None,
    group_id: str | int | None = None,
) -> list[Task]:
    """Set up a group of subtasks from a run dictionary, creating run directories and
    returning them in a table

    """
    subpes_dct = util.subpes_dict_from_mechanism_dat(file_dct.get("mechanism.dat"))
    keys = util.subtask_keys_from_run_dict(
        run_dct, task_type, key_type, subpes_dct=subpes_dct
    )
    task_lines = util.task_lines_from_run_dict(run_dct, task_type, key_type)
    tasks = []
    for task_key, task_line in enumerate(task_lines):
        task_name = util.parse_task_name(task_line)
        task_path = Path(f"{group_id}_{task_key:02d}_{task_name}")
        nworkers_lst = util.parse_subtasks_nworkers(
            task_line, file_dct, keys, task_type=task_type, key_type=key_type
        )
        if "pes_groups.dat" in file_dct and task_type == "ktp":
            subtasks = [
                Subtask(
                    key=run_dct.get("pes"),
                    nworkers=nworkers_lst[0],
                    path=str(task_path / "all"),
                )
            ]
        else:
            subtasks = [
                Subtask(
                    key=key,
                    nworkers=nworkers,
                    path=str(task_path / util.format_subtask_key(key)),
                )
                for key, nworkers in zip(keys, nworkers_lst, strict=True)
            ]
        tasks.append(
            Task(
                name=task_name,
                line=task_line,
                mem=util.parse_task_memory(task_line, file_dct),
                nprocs=util.parse_task_nprocs(task_line, file_dct),
                subtasks=subtasks,
            )
        )

    if tasks and task_type in COMBINED_TASK_GROUPS:
        idx = next((i for i, t in enumerate(tasks) if t.name == "run_mess"), 0)
        task = tasks[idx]
        task.line = "\n".join(t.line for t in tasks)
        tasks = [task]

    return tasks
