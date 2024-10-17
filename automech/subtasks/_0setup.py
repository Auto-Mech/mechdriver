""" Standalone script to break an AutoMech input into subtasks for parallel execution
"""

import dataclasses
import io
import os
import re
import textwrap
from collections.abc import Sequence
from pathlib import Path

import automol
import more_itertools as mit
import pandas
import pyparsing as pp
import yaml
from pyparsing import common as ppc

COMMENT_REGEX = re.compile(r"#.*$", flags=re.M)
ALL_KEY = "all"

DEFAULT_MEM = 20
DEFAULT_TASK_GROUPS = ("els", "thermo", "ktp")
GROUP_ID = {"els-spc": 0, "els-pes": 1, "thermo": 2, "ktp": 3}
GROUP_TASK_AND_KEY_TYPE = {
    "els-spc": ("els", "spc"),
    "els-pes": ("els", "pes"),
    "thermo": ("thermo", "spc"),
    "ktp": ("ktp", None),
}

SUBTASK_DIR = "subtasks"
INFO_FILE = "info.yaml"

ROTOR_TASKS = ("hr_scan",)
SAMP_TASKS = ("conf_samp",)


class TableKey:
    task = "task"


class InfoKey:
    group_ids = "group_ids"  # identifiers for each subtask group, in the order they should be run
    run_path = "run_path"  # path to the run filesystem
    save_path = "save_path"  # path to the save filesystem
    # If relative paths are given for `path`, `save_path`, and/or `run_path`, they will
    # be relative to the following `work_path`, which is the user's current working
    # directory when they run the setup command:
    work_path = "work_path"  # path to where the user ran the command


@dataclasses.dataclass
class Task:
    name: str
    line: str
    mem: int
    nprocs: int
    subtask_keys: list[str]
    subtask_nworkers: list[int]


def setup(
    path: str | Path,
    out_path: str | Path = SUBTASK_DIR,
    save_path: str | Path | None = None,
    run_path: str | Path | None = None,
    task_groups: Sequence[str] = DEFAULT_TASK_GROUPS,
):
    """Creates run directories for each task/species/TS and returns the paths in tables

    Task types: 'els', 'thermo', or 'ktp'
    Subtask types: 'spc', 'pes', or `None` (=all species and/or reactions)

    :param path: The path to the AutoMech input to split into subtasks
    :param out_path: The root path of the output (will be filled with subtask
        directories, CSVs, and YAML file)
    :param save_path: The path to the save filesystem
        (if `None`, the value in run.dat is used)
    :param run_path: The path to the run filesystem
        (if `None`, the value in run.dat is used)
    :param task_groups: The task groups to set up
    :return: DataFrames of run paths, whose columns (species/TS index) are independent
        and can be run in parallel, but whose rows (tasks) are potentially sequential
    """
    # Pre-process the task groups
    task_groups = list(task_groups)
    if "els" in task_groups:
        idx = task_groups.index("els")
        task_groups[idx : idx + 1] = ("els-spc", "els-pes")
    assert all(g in GROUP_ID for g in task_groups), f"{task_groups} not in {GROUP_ID}"

    # Read input files from source path
    path = Path(path)
    file_dct = read_input_files(path)
    run_dct = parse_run_dat(file_dct.get("run.dat"))

    # Set the run and save paths
    save_path, run_path = filesystem_paths_from_run_dict(
        run_dct, save_path=save_path, run_path=run_path
    )
    run_dct["input"] = f"run_prefix = {run_path}\nsave_prefix = {save_path}"

    # Create the path for the subtask directories
    out_path = Path(out_path)
    out_path.mkdir(exist_ok=True)

    # Set up the subtasks for each group
    run_group_ids = []
    for task_group in task_groups:
        group_id = GROUP_ID.get(task_group)
        task_type, key_type = GROUP_TASK_AND_KEY_TYPE.get(task_group)
        ret = setup_subtask_group(
            run_dct,
            file_dct,
            task_type=task_type,
            key_type=key_type,
            group_id=group_id,
            out_path=out_path,
        )
        if ret is not None:
            run_group_ids.append(group_id)

    # Write the subtask info to YAML
    info_path = out_path / INFO_FILE
    print(f"Writing general information to {info_path}")
    info_dct = {
        InfoKey.save_path: save_path,
        InfoKey.run_path: run_path,
        InfoKey.work_path: os.getcwd(),
        InfoKey.group_ids: run_group_ids,
    }
    info_path.write_text(yaml.safe_dump(info_dct))


def setup_subtask_group(
    run_dct: dict[str, str],
    file_dct: dict[str, str],
    task_type: str,
    key_type: str | None = None,
    group_id: str | int | None = None,
    out_path: str | Path = SUBTASK_DIR,
) -> pandas.DataFrame | None:
    """Set up a group of subtasks from a run dictionary, creating run directories and
    returning them in a table

    :param source_path: The path to the AutoMech input to split into subtasks
    :param task_type: The type of task: 'els', 'thermo', or 'ktp'
    :param key_type: The type of subtask key: 'spc', 'pes', or `None`
    :param group_id: The group ID, used to name files and folders
    :return: A DataFrame of run paths, whose columns (subtasks) are independent and can
        be run in parallel, but whose rows (tasks) are potentially sequential
    """

    def _subtask_directory(key):
        id_str_ = "{:02d}".format
        key_ = parse_subtask_key(key)
        subtask_dir = key_
        if isinstance(key_, tuple):
            subtask_dir = "_".join(map(id_str_, key_))
        if isinstance(key_, int):
            subtask_dir = id_str_(key_)
        assert isinstance(subtask_dir, str), f"Invalid subtask key: {key}"
        return subtask_dir

    # Form a prefix for the task/subtask type
    if group_id is None:
        type_keys = [task_type] + ([] if key_type is None else [key_type])
        return "_".join(type_keys)
    group_id = str(group_id)

    # Blocks that must be included in the run.dat
    block_keys = ["input"] + (["pes", "spc"] if key_type is None else [key_type])

    # Parse tasks and subtask keys for this group
    tasks = determine_task_list(run_dct, file_dct, task_type, key_type)

    # If the task list is empty, return `None`
    if not tasks:
        return None

    # Get the specs for each task and write them to a YAML file
    yaml_path = out_path / f"{group_id}.yaml"
    print(f"Writing task specs to {yaml_path}")
    write_task_list(tasks, yaml_path)

    # Create directories for each subtask and save the paths in a DataFrame
    row_dcts = []
    for task_key, task in enumerate(tasks):
        row_dct = {TableKey.task: task.name}
        task_path = out_path / f"{group_id}_{task_key:02d}_{task.name}"
        print(f"Setting up subtask directories in {task_path}")
        task_run_dct = {k: v for k, v in run_dct.items() if k in block_keys}
        task_run_dct[task_type] = task.line
        for key in task.subtask_keys:
            # Generate the subtask path
            subtask_path = task_path / _subtask_directory(key)
            subtask_path.mkdir(parents=True, exist_ok=True)
            # Generate the input file dictionary
            subtask_run_dct = task_run_dct.copy()
            if key != ALL_KEY:
                subtask_run_dct[key_type] = key
            subtask_file_dct = {**file_dct, "run.dat": form_run_dat(subtask_run_dct)}
            # Write the input files and append the path to the current dataframe row
            write_input_files(subtask_path, subtask_file_dct)
            row_dct[key] = subtask_path
        row_dcts.append(row_dct)

    df = pandas.DataFrame(row_dcts)

    # Write the subtask table to CSV
    csv_path = out_path / f"{group_id}.csv"
    print(f"Writing subtask table to {csv_path}")
    df.to_csv(csv_path, index=False)
    print()

    return df


def determine_task_list(
    run_dct: dict[str, str],
    file_dct: dict[str, str],
    task_type: str,
    key_type: str | None = None,
) -> list[Task]:
    """Set up a group of subtasks from a run dictionary, creating run directories and
    returning them in a table

    """
    keys = subtask_keys_from_run_dict(run_dct, key_type)

    tasks: list[Task] = [
        Task(
            name=parse_task_name(task_line),
            line=task_line,
            mem=parse_task_memory(task_line, file_dct),
            nprocs=parse_task_nprocs(task_line, file_dct),
            subtask_keys=keys,
            subtask_nworkers=parse_subtasks_nworkers(
                task_line, file_dct, subtask_keys=keys
            ),
        )
        for task_line in task_lines_from_run_dict(run_dct, task_type, key_type)
    ]

    fit_idx = next((i for i, t in enumerate(tasks) if t.name == "run_fits"), None)
    run_idx = next((i for i, t in enumerate(tasks) if t.name == "run_mess"), None)

    if fit_idx is not None and run_idx is not None:
        fit_task = tasks.pop(fit_idx)
        tasks[run_idx].line += f"\n{fit_task.line}"

    return tasks


# Functions acting on the run directory as a whole
def read_input_files(run_dir: str | Path) -> dict[str, str]:
    inp_dir = Path(run_dir) / "inp"
    return {
        "run.dat": (inp_dir / "run.dat").read_text(),
        "theory.dat": (inp_dir / "theory.dat").read_text(),
        "models.dat": (inp_dir / "models.dat").read_text(),
        "mechanism.dat": (inp_dir / "mechanism.dat").read_text(),
        "species.csv": (inp_dir / "species.csv").read_text(),
    }


def write_input_files(run_dir: str | Path, file_dct: dict[str, str]) -> None:
    inp_dir = Path(run_dir) / "inp"
    inp_dir.mkdir(exist_ok=True)
    for name, contents in file_dct.items():
        (inp_dir / name).write_text(contents)


# Parse task information
def parse_task_name(task_line: str) -> str:
    """Parse the task name from a task line

    :param task_line: The task line from the run.dat file
    :return: The task name
    """
    task_name = task_line.split()[0]
    if task_name in ["spc", "ts", "all"]:
        task_name = task_line.split()[1]
    return task_name


def parse_task_fields(task_line: str) -> dict[str, str]:
    """Parse in-line fields of the form `key=val`

    :param inp: The string to parse
    :return: The fields, as a dictionary
    """
    word = pp.Word(pp.printables, exclude_chars="=")
    eq = pp.Suppress(pp.Literal("="))
    field = pp.Group(word + eq + word)
    expr = pp.Suppress(...) + pp.DelimitedList(field, delim=pp.WordEnd()) | pp.Empty()
    field_dct: dict[str, str] = dict(expr.parseString(task_line).as_list())
    return field_dct


def parse_task_memory(task_line: str, file_dct: dict[str, str]) -> int:
    """Parse the memory requirement for a given task

    :param task_line: The task line from the run.dat file
    :param file_dct: The file dictionary
    :return: The memory requirement for the task
    """
    field_dct = parse_task_fields(task_line)

    mem = DEFAULT_MEM

    if "runlvl" in field_dct:
        runlvl = field_dct.get("runlvl")
        theory_dct = parse_theory_dat(file_dct.get("theory.dat"))
        mem = int(float(theory_dct.get(runlvl).get("mem")))

    return mem


def parse_task_nprocs(task_line: str, file_dct: dict[str, str]) -> int:
    """Read the memory and nprocs specs for a given task

    :param task_line: The task line from the run.dat file
    :param file_dct: The file dictionary
    :return: The memory and nprocs for the task
    """
    field_dct = parse_task_fields(task_line)

    nprocs = 1

    if "nprocs" in field_dct:
        nprocs = int(float(field_dct.get("nprocs")))

    if "runlvl" in field_dct:
        runlvl = field_dct.get("runlvl")
        theory_dct = parse_theory_dat(file_dct.get("theory.dat"))
        nprocs = int(float(theory_dct.get(runlvl).get("nprocs")))

    return nprocs


def parse_subtasks_nworkers(
    task_line: str, file_dct: dict[str, str], subtask_keys: list[str]
) -> list[int]:
    """Read the memory and nprocs specs for a given task

    :param task_line: The task line from the run.dat file
    :param file_dct: The file dictionary
    :param nsub: The
    :return: The memory and nprocs for the task
    """
    nworkers_lst = [1] * len(subtask_keys)

    if task_line.startswith("spc"):
        task_name = parse_task_name(task_line)
        field_dct = parse_task_fields(task_line)

        # Get the list of ChIs ordered by subtask key
        spc_df = parse_species_csv(file_dct.get("species.csv"))
        if "inchi" not in spc_df:
            spc_df["inchi"] = spc_df["smiles"].apply(automol.smiles.inchi)
        chis = [spc_df.iloc[int(k) - 1]["inchi"] for k in subtask_keys]

        # Determine the number of workers per subtask
        if task_name in ROTOR_TASKS:
            nworkers_lst = list(map(rotor_count_from_inchi, chis))
        if task_name in SAMP_TASKS or field_dct.get("cnf_range", "").startswith("n"):
            nmax = int(field_dct.get("cnf_range", "n100")[1:])
            nsamp_lst = [sample_count_from_inchi(c, param_d=nmax) for c in chis]
            nworkers_lst = [max((n - 1) // 2, 1) for n in nsamp_lst]

    return nworkers_lst


# Functions acting on theory.dat data
def parse_theory_dat(theory_dat: str) -> dict[str, dict[str, str]]:
    """Parse a theory.dat file into a dictionary of dictionaries

    :param theory_dat: The contents of the theory.dat file, as a string
    :return: The dictionary of the parsed theory.dat file
    """
    theory_dat = without_comments(theory_dat)
    theory_expr = pp.OneOrMore(block_expression("level", key="content"))
    blocks = theory_expr.parseString(theory_dat).as_list()

    word = pp.Word(pp.printables)
    eq = pp.Suppress(pp.Literal("="))
    field = pp.Group(word("key") + eq + word("val"))
    block_expr = word("key") + pp.DelimitedList(field, delim=pp.LineEnd())("fields")

    theory_dct = {}
    for block in blocks:
        res = block_expr.parseString(block)
        theory_dct[res.get("key")] = dict(res.get("fields").as_list())

    return theory_dct


# Functions acting on species.csv data
def parse_species_csv(species_csv: str) -> pandas.DataFrame:
    """Parse a species.csv file into a pandas dataframe

    :param species_csv: The contents of the species.csv file, as a string
    :return: The species table
    """
    return pandas.read_csv(io.StringIO(species_csv), quotechar="'")


def rotor_count_from_inchi(chi: str) -> int:
    """Determine species rotor count for a molecule from its InChI or AMChI string

    :param chi: An InChI or AMChI string
    :return: The rotor count
    """
    gra = automol.amchi.graph(chi, stereo=False)
    # If there are no torsions at all, return 1
    if not len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=True)):
        return 1

    nrotor = len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=True))
    return nrotor


def sample_count_from_inchi(
    chi: str,
    param_a: int = 12,
    param_b: int = 1,
    param_c: int = 3,
    param_d: int = 100,
) -> int:
    """Determine species MC sample count for a molecule from its InChI or AMChI string

    The parameters (a, b, c, d) are used to calculate the sample count as follows:
    ```
        nsamp = min(a + b * c^ntors, d)
    ```
    where `ntors` is the number of torsional degrees of freedom in the molecule.

    :param chi: An InChI or AMChI string
    :param param_a: The `a` parameter used to calculate the sample count
    :param param_b: The `b` parameter used to calculate the sample count
    :param param_c: The `c` parameter used to calculate the sample count
    :param param_d: The `d` parameter used to calculate the sample count
    :return: The sample count
    """
    gra = automol.amchi.graph(chi, stereo=False)
    # If there are no torsions at all, return 1
    if not len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=True)):
        return 1

    ntors = len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=False))
    nsamp = min(param_a + param_b * param_c**ntors, param_d)
    return nsamp


# Functions acting on run.dat data
def parse_run_dat(run_dat: str) -> dict[str, str]:
    """Parse a run.dat file into a dictionary of blocks

    :param run_dat: The contents of the run.dat file, as a string
    :return: The dictionary of the parsed run.dat file
    """

    def _parse_block(run_dat, keyword):
        expr = block_expression(keyword, key="content")
        res, *_ = next(expr.scan_string(run_dat), [None])
        if res is None:
            return None

        content = res.get("content")
        return format_block(content)

    run_dat = without_comments(run_dat)
    block_dct = {
        "input": _parse_block(run_dat, "input"),
        "pes": _parse_block(run_dat, "pes"),
        "spc": _parse_block(run_dat, "spc"),
        "els": _parse_block(run_dat, "els"),
        "thermo": _parse_block(run_dat, "thermo"),
        "ktp": _parse_block(run_dat, "ktp"),
    }
    return {k: v for k, v in block_dct.items() if v is not None}


def form_run_dat(run_dct: dict[str, str]) -> str:
    """Format the contents of a run.dat file from a run.dat dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :return: The run.dat file contents, as a string
    """
    keys = ["input", "spc", "pes", "els", "thermo", "ktp"]
    run_dat = ""
    for key in keys:
        if key in run_dct:
            run_dat += f"{key}\n{format_block(run_dct.get(key))}\nend {key}\n\n"
    return run_dat


def filesystem_paths_from_run_dict(
    run_dct: dict[str, str],
    save_path: str | Path | None = None,
    run_path: str | Path | None = None,
) -> str:
    """Get the input block of a run dictionary, with absolute paths for the RUN and SAVE
    directories

    :param run_dct: The dictionary of a parsed run.dat file
    :param save_path: The path to the save filesystem
        (if `None`, the value in run.dat is used.)
    :param run_path: The path to the run filesystem
        (if `None`, the value in run.dat is used.)
    :return: The input block, with absolute paths
    """

    def _extract_path(key: str) -> str:
        inp_block = run_dct.get("input")
        inp_block = without_comments(inp_block)
        word = pp.Word(pp.printables, exclude_chars="=")
        field = pp.Suppress(... + pp.Literal(f"{key}_prefix") + pp.Literal("="))
        expr = field + word("path")
        return expr.parseString(inp_block).get("path")

    save_path = _extract_path("save") if save_path is None else save_path
    run_path = _extract_path("run") if run_path is None else run_path
    return save_path, run_path


def subtask_keys_from_run_dict(
    run_dct: dict[str, str], subtask_type: str | None = None
) -> list[str]:
    """Extract species indices from a run.dat dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :return: A sequence of indices for each individual species
    """

    if subtask_type is None:
        return [ALL_KEY]

    if subtask_type == "spc" and "spc" in run_dct:
        spc_block = run_dct.get("spc")
        return list(map(str, parse_index_series(spc_block)))

    if subtask_type == "pes" and "pes" in run_dct:
        pes_block = run_dct.get("pes")

        colon = pp.Suppress(pp.Literal(":"))
        before_line_end = pp.FollowedBy(pp.LineEnd())
        entry = pp.Group(
            ppc.integer("pes") + colon + pp.SkipTo(before_line_end)("channels")
        )
        expr = pp.DelimitedList(entry, delim=pp.LineEnd())

        keys = []
        for res in expr.parseString(pes_block):
            pidx = res.get("pes")
            cidxs = parse_index_series(res.get("channels"))
            keys.extend(f"{pidx}: {cidx}" for cidx in cidxs)
        return list(mit.unique_everseen(keys))

    return []


def task_lines_from_run_dict(
    run_dct: dict[str, str], task_type: str, subtask_type: str | None = None
) -> list[str]:
    """Extract electronic structure tasks from  of a run.dat dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :param task_type: The type of task: 'els', 'thermo', or 'ktp'
    :param subtask_type: The type of subtask: 'spc', 'pes', or `None`
    :return: A sequence of (task name, task line) pairs
    """
    if task_type not in run_dct:
        return []

    block = run_dct.get(task_type)
    lines = [line.strip() for line in block.splitlines()]
    if task_type == "els":
        types = ("spc", "pes")
        assert subtask_type in types, f"Subtask type {subtask_type} not in {types}"
        start_key = "ts" if subtask_type == "pes" else "spc"
        lines = [
            line
            for line in lines
            if line.startswith(start_key) or line.startswith("all")
        ]

    return lines


# Task read/write
def write_task_list(yaml_tasks: list[Task], path: str | Path) -> None:
    """Write a task list out in YAML format

    :param tasks: The list of tasks
    :param path: The path to the YAML file to write
    """
    path = Path(path)
    yaml_tasks = list(map(dataclasses.asdict, yaml_tasks))
    path.write_text(yaml.safe_dump(yaml_tasks, default_flow_style=None))


def read_task_list(path: str | Path) -> list[Task]:
    """Read a list of tasks from a YAML file

    :param path: The path to the YAML file
    :return: The list of tasks
    """
    path = Path(path)
    yaml_tasks = yaml.safe_load(path.read_text())
    return [Task(**d) for d in yaml_tasks]


# Generic string formatting functions
def format_block(inp: str) -> str:
    """Format a block with nice indentation

    :param inp: A multiline string to be formatted
    :return: The formatted string
    """
    inp = re.sub(r"^\ *", "", inp, flags=re.MULTILINE)
    return textwrap.indent(inp, "    ")


def without_comments(inp: str) -> str:
    """Get a CHEMKIN string or substring with comments removed.

    :param inp: A CHEMKIN mechanism, as a file path or string
    :return: The string, without comments
    """
    return re.sub(COMMENT_REGEX, "", inp)


def parse_subtask_key(key: str) -> int | tuple[int, int] | str:
    """Parse a subtask key into its components

    Examples:
        '1'     ->  1
        '1: 2'  ->  (1, 2)
        'all'   ->  'all'

    :param key: The key to parse
    :return: The parsed components
    """
    spc_key = ppc.integer
    all_key = pp.Literal(ALL_KEY)
    pes_key = ppc.integer + pp.Suppress(":") + ppc.integer
    expr = (spc_key ^ all_key ^ pes_key) + pp.StringEnd()
    res = expr.parseString(key).as_list()
    return res[0] if len(res) == 1 else tuple(res)


def parse_index_series(inp: str) -> list[int]:
    r"""Parse a sequence of indices from a string separated by commas and newlines,
    with ranges indicated by 'x-y'

    Example:
        Input: '1,3, 5-9  \n  11,13-14\n23\n 27-29'
        Output: (1, 3, 5, 6, 7, 8, 9, 11, 13, 14, 23, 27, 28, 29)
    """
    if not inp:
        return []

    dash = pp.Suppress(pp.Literal("-"))
    entry = ppc.integer ^ pp.Group(ppc.integer + dash + ppc.integer)
    delim = pp.LineEnd() ^ pp.Literal(",")
    expr = pp.DelimitedList(entry, delim=delim)
    idxs = []
    for res in expr.parseString(inp).as_list():
        if isinstance(res, int):
            idxs.append(res)
        else:
            start, stop = res
            idxs.extend(range(start, stop + 1))
    return list(mit.unique_everseen(idxs))


def block_expression(keyword: str, key: str = "content") -> pp.ParseExpression:
    """Parse a block from an AutoMech input file

    :param inp: The input file contents
    :param keyword: The block keyword
    :return: _description_
    """
    start = pp.Keyword(keyword)
    end = pp.Keyword("end") + pp.Keyword(keyword)
    return pp.Suppress(... + start) + pp.SkipTo(end)(key) + pp.Suppress(end)
