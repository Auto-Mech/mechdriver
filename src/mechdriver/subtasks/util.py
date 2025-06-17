"""Standalone script to break an AutoMech input into subtasks for parallel execution"""

import io
import itertools
import re
import textwrap
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

import automol
import more_itertools as mit
import pandas
import pyparsing as pp
from pyparsing import common as ppc

COMMENT_REGEX = re.compile(r"#.*$", flags=re.M)
ALL_KEY = "all"

DEFAULT_MEM = 20


ROTOR_TASKS = ("hr_scan",)
SAMP_TASKS = ("conf_samp",)
INP_DIR_TASKS = ("conf_energy",)


# Functions acting on the run directory as a whole
def read_input_files(run_dir: str | Path) -> dict[str, str]:
    inp_dir = Path(run_dir) / "inp"
    return {
        fname: (inp_dir / fname).read_text() 
        for fname in [
            "run.dat", "theory.dat", "models.dat", 
            "mechanism.dat", "species.csv", "pes_groups.dat"
        ] 
        if (inp_dir / fname).exists()
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
        lvl = field_dct.get("runlvl")
        theory_dct = parse_theory_dat(file_dct.get("theory.dat"))
        mem = int(float(theory_dct.get(lvl).get("mem")))

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
        lvl = field_dct.get("runlvl")
        theory_dct = parse_theory_dat(file_dct.get("theory.dat"))
        nprocs = int(float(theory_dct.get(lvl).get("nprocs")))

    return nprocs


def parse_subtasks_nworkers(
    task_line: str,
    file_dct: dict[str, str],
    subtask_keys: list[str],
    task_type: str,
    key_type: str | None = None,
) -> list[int]:
    """Read the memory and nprocs specs for a given task

    :param task_line: The task line from the run.dat file
    :param file_dct: The file dictionary
    :param nsub: The
    :return: The memory and nprocs for the task
    """
    nworkers_lst = [1] * len(subtask_keys)

    if task_type != "els":
        return nworkers_lst

    task_name = parse_task_name(task_line)
    field_dct = parse_task_fields(task_line)

    # Get the list of ChIs ordered by subtask key
    spc_df = parse_species_csv(file_dct.get("species.csv"))

    gras_lst = None
    if key_type == "spc":
        # Get the list of ChIs ordered by subtask key
        chis = [spc_df.iloc[int(k) - 1]["inchi"] for k in subtask_keys]
        gras_lst = [[automol.amchi.graph(c, stereo=False)] for c in chis]

    if key_type == "pes":
        task_name = parse_task_name(task_line)
        field_dct = parse_task_fields(task_line)

        chi_dct = dict(zip(spc_df["name"], spc_df["inchi"], strict=True))
        rxn_dct = parse_mechanism_dat(file_dct.get("mechanism.dat"))
        if rxn_dct:
            rxns = [rxn_dct.get(k) for k in subtask_keys]
            chis = [tuple(list(map(chi_dct.get, rgts)) for rgts in rxn) for rxn in rxns]
            gras_lst = [
                list(
                    map(automol.reac.ts_graph, automol.reac.from_chis(*c, stereo=False))
                )
                for c in chis
            ]

    if gras_lst is None:
        return nworkers_lst

    # Determine the number of workers per subtask
    if task_name in ROTOR_TASKS:
        nworkers_lst = list(map(rotor_count_from_graphs, gras_lst))
    if task_name in SAMP_TASKS or field_dct.get("cnf_range", "").startswith("n"):
        nmax = int(field_dct.get("cnf_range", "n100").split('_')[0][1:])
        nsamp_lst = [sample_count_from_graphs(gs, param_d=nmax) for gs in gras_lst]
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
    spc_df = pandas.read_csv(io.StringIO(species_csv), quotechar="'")

    if "inchi" not in spc_df:
        spc_df["inchi"] = spc_df["smiles"].apply(automol.smiles.chi)

    #if "canon_enant_ich" not in spc_df:
    #    for spc in spc_df.index:
    #        try:
    #            spc_df.loc[spc, "canon_enant_ich"] = automol.chi.canonical_enantiomer(spc_df["inchi"][spc])
    #        except AssertionError:
    #           spc_df.loc[spc, "canon_enant_ich"] = spc_df["inchi"][spc]

    return spc_df


def rotor_count_from_graphs(gras: Sequence[object]) -> int:
    """Determine the total rotor count for a sequence of graphs

    :param gras: A sequence of molecule or TS graphs
    :return: The rotor count
    """
    return sum(map(rotor_count_from_graph, gras))


def sample_count_from_graphs(
    gras: Sequence[object],
    param_a: int = 12,
    param_b: int = 1,
    param_c: int = 3,
    param_d: int = 100,
) -> int:
    """Determine the total conformer sample count for a sequence of graphs

    The parameters (a, b, c, d) are used to calculate the sample count as follows:
    ```
        nsamp = min(a + b * c^ntors, d)
    ```
    where `ntors` is the number of torsional degrees of freedom in the molecule.

    :param gra: A molecule or TS graph
    :param param_a: The `a` parameter used to calculate the sample count
    :param param_b: The `b` parameter used to calculate the sample count
    :param param_c: The `c` parameter used to calculate the sample count
    :param param_d: The `d` parameter used to calculate the sample count
    :return: The sample count
    """
    return sum(
        sample_count_from_graph(
            g, param_a=param_a, param_b=param_b, param_c=param_c, param_d=param_d
        )
        for g in gras
    )


def rotor_count_from_graph(gra: object) -> int:
    """Determine the rotor count for a graph

    :param gra: A molecule or TS graph
    :return: The rotor count
    """
    # If there are no torsions at all, return 1
    if not len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=True)):
        return 1

    nrotor = len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=True))
    return nrotor


def sample_count_from_graph(
    gra: object,
    param_a: int = 12,
    param_b: int = 1,
    param_c: int = 3,
    param_d: int = 100,
) -> int:
    """Determine the conformer sample count for a graph

    The parameters (a, b, c, d) are used to calculate the sample count as follows:
    ```
        nsamp = min(a + b * c^ntors, d)
    ```
    where `ntors` is the number of torsional degrees of freedom in the molecule.

    :param gra: A molecule or TS graph
    :param param_a: The `a` parameter used to calculate the sample count
    :param param_b: The `b` parameter used to calculate the sample count
    :param param_c: The `c` parameter used to calculate the sample count
    :param param_d: The `d` parameter used to calculate the sample count
    :return: The sample count
    """
    # If there are no torsions at all, return 1
    if not len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=True)):
        return 1

    ntors = len(automol.graph.rotational_bond_keys(gra, with_ch_rotors=False))
    nsamp = min(param_a + param_b * param_c**ntors, param_d)
    return nsamp


# Functions acting on other chemical data
def ts_multiplicity(
    rct_muls: Sequence[int], prd_muls: Sequence[int], high: bool = False
) -> int:
    """Guess the TS multiplicity from reactants and products.

    :param rct_muls: Reactant multiplicities
    :param prd_muls: Product multiplicities
    :param high: Whether to guess the high-spin multiplicity
    :return: TS multiplicity
    """
    if len(rct_muls) == len(prd_muls) == 1:
        return max(rct_muls[0], prd_muls[0])

    rct_spin_sum = sum(m - 1 for m in rct_muls)
    prd_spin_sum = sum(m - 1 for m in prd_muls)
    ts_spin = (
        max(rct_spin_sum, prd_spin_sum) if high else min(rct_spin_sum, prd_spin_sum)
    )
    ts_mul = int(round(ts_spin + 1))
    return ts_mul


# Functions acting on mechanism.dat data
def parse_mechanism_dat(mechanism_dat: str) -> dict[str, tuple[list[str], list[str]]]:
    """Parse the mechanism.dat file into a list of dictionaries.

    {
        "1: 2": (["C2H6", "OH"], ["C2H5", "H2O"]),
        ...
    }

    :param mechanism_dat: The contents of the mechanism.dat file, as a string
    :return: A list of information for each reaction
    """

    class Key:
        eq = "eq"
        arrh = "arrh"
        comment = "comment"

    sort_key = sort_val = pp.DelimitedList(
        pp.Word(pp.alphanums + "_", exclude_chars="."), delim=".", min=3
    )
    sort_expr = pp.Group(sort_key) + pp.Group(sort_val)

    comm_mark = pp.Suppress(pp.Char("!") ^ pp.Char("#"))
    comm_expr = pp.SkipTo(pp.LineEnd())
    arrh_expr = pp.WordStart() + ppc.number[3]
    eq_expr = pp.SkipTo(arrh_expr).set_parse_action(lambda t: str.rstrip(t[0]))
    reac_expr = (
        eq_expr(Key.eq) + arrh_expr(Key.arrh) + comm_mark + comm_expr(Key.comment)
    )

    reac_block = re.search(
        "REACTIONS(.*?)END", mechanism_dat, flags=re.M | re.S | re.I
    ).group(1)
    reac_dct = {}
    for line in reac_block.splitlines():
        if reac_expr.matches(line):
            res = reac_expr.parse_string(line)
            eq = res.get(Key.eq)
            comment = res.get(Key.comment)
            sort_info = dict(
                itertools.zip_longest(*sort_expr.parse_string(comment).as_list(),
                                      fillvalue='MISSING')
            )
            pes = int(sort_info.get("pes"))
            channel = int(sort_info.get("channel"))
            reac_dct[f"{pes}: {channel}"] = tuple(
                [s.strip() for s in re.split(r"\+(?!\s*\+)", r)]
                for r in re.split("=|=>|<=>", eq)
            )

    return reac_dct


def subpes_dict_from_mechanism_dat(
    mechanism_dat: str,
) -> dict[str, list[list[str]]] | None:
    """Determine the groups of channels for each sub-PES from the mechanism.dat file

    Structure of the return dictionary:

        {'1': [['1', '2', '3'], ['4', '5', '6', '7', '8']],
         '2': [['1'], ['2'], ['3']],
         '3': [['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']],
         ...}

    where the keys are PES indices and the values are grouped channel indices.

    :param mechanism_dat: The contents of the mechanism.dat file, as a string
    :return: A dictionary mapping PESs onto groups of channels, by string index
    """
    # Define a parser to extract "pes.subpes.channel 1.2.3" data
    comment_mark = pp.Char("!") | pp.Char("#")
    sort_key = sort_val = pp.DelimitedList(
        pp.Word(pp.alphanums + "_"), delim=".", combine=True, min=3
    )
    sort_item = pp.Suppress(...) + pp.Group(
        pp.Suppress(comment_mark) + sort_key + sort_val + pp.Suppress(pp.LineEnd())
    )
    sort_items = pp.ZeroOrMore(sort_item)
    results: list[str, str] = sort_items.parse_string(mechanism_dat).as_list()
    if not results:
        return None

    # Create a DataFrame with "pes | subpes | channel" columns
    sort_df = pandas.DataFrame.from_records(
        [dict(zip(k.split("."), v.split("."), strict=True)) for k, v in results]
    )
    sort_df = sort_df.apply(pandas.to_numeric, axis=1, errors='coerce')

    if not all(k in sort_df for k in ("pes", "subpes", "channel")):
        return None

    # Create a dictionary with groups of channels by subpes
    subpes_dct = defaultdict(list)
    for (p, _), d in sort_df.groupby(["pes", "subpes"])["channel"]:
        subpes_dct[p].append(d.to_list())

    return dict(subpes_dct)


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
        "proc": _parse_block(run_dat, "proc"),
    }

    return {k: v for k, v in block_dct.items() if v is not None}


def form_run_dat(run_dct: dict[str, str]) -> str:
    """Format the contents of a run.dat file from a run.dat dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :return: The run.dat file contents, as a string
    """
    keys = ["input", "spc", "pes", "els", "thermo", "ktp", "proc"]
    run_dat = ""
    for key in keys:
        if key in run_dct:
            run_dat += f"{key}\n{format_block(run_dct.get(key))}\nend {key}\n\n"

    return run_dat


def input_arguments_from_run_dict(
    run_dct: dict[str, str],
    save_path: str | Path | None = None,
    run_path: str | Path | None = None,
    absolute: bool = True,
) -> dict[str, str]:
    """Get the input block of a run dictionary, with absolute paths for the RUN and SAVE
    directories.

    :param run_dct: The dictionary of a parsed run.dat file
    :param save_path: The path to the save filesystem
        (if `None`, the value in run.dat is used.)
    :param run_path: The path to the run filesystem
        (if `None`, the value in run.dat is used.)
    :param absolute: Whether to resolve absolute paths
    :return: The input block, with absolute paths
    """
    inp_block = run_dct.get("input")
    inp_block = without_comments(inp_block)
    lines = inp_block.strip().splitlines()
    inp_dct = dict(
        line.replace(" ", "").split("=", maxsplit=1) for line in lines if "=" in line
    )
    save_path = save_path or inp_dct.get("save_prefix")
    run_path = run_path or inp_dct.get("run_prefix")
    if absolute:
        save_path = str(Path(save_path).resolve())
        run_path = str(Path(run_path).resolve())

    inp_dct["save_prefix"] = save_path
    inp_dct["run_prefix"] = run_path
    return inp_dct


def filesystem_paths_from_run_dict(
    run_dct: dict[str, str],
    save_path: str | Path | None = None,
    run_path: str | Path | None = None,
) -> tuple[str, str]:
    """Get the input block of a run dictionary, with absolute paths for the RUN and SAVE
    directories.

    :param run_dct: The dictionary of a parsed run.dat file
    :param save_path: The path to the save filesystem
        (if `None`, the value in run.dat is used.)
    :param run_path: The path to the run filesystem
        (if `None`, the value in run.dat is used.)
    :return: The input block, with absolute paths
    """
    inp_dct = input_arguments_from_run_dict(run_dct, save_path=save_path, run_path=run_path)
    return inp_dct.get("save_prefix"), inp_dct.get("run_prefix")


def subtask_keys_from_run_dict(
    run_dct: dict[str, str],
    task_type: str,
    key_type: str | None = None,
    subpes_dct: dict[str, list[list[str]]] | None = None,
) -> list[str]:
    """Extract species indices from a run.dat dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :return: A sequence of indices for each individual species
    """

    if key_type is None:
        return [ALL_KEY]

    if key_type == "spc" and "spc" in run_dct and run_dct.get("spc").strip():
        spc_block = run_dct.get("spc")
        return list(map(str, parse_index_series(spc_block)))

    if key_type == "pes" and "pes" in run_dct and run_dct.get("pes").strip():
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
            cidx_range = res.get("channels")
            cidxs = parse_index_series(cidx_range)
            if task_type == "ktp" and subpes_dct is None:
                keys.append(f"{pidx}: {cidx_range}")
            elif task_type == "ktp":
                assert pidx in subpes_dct
                cidx_groups = [
                    (x for x in xs if x in cidxs) for xs in subpes_dct.get(pidx)
                ]
                keys.extend(
                    f"{pidx}: {format_index_series(cidx_group)}"
                    for cidx_group in cidx_groups
                )
            else:
                keys.extend(f"{pidx}: {cidx}" for cidx in cidxs)
        return list(mit.unique_everseen(keys))

    return []


def task_lines_from_run_dict(
    run_dct: dict[str, str], task_type: str, key_type: str | None = None
) -> list[str]:
    """Extract tasks from a run.dat dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :param task_type: The type of task: 'els', 'thermo', or 'ktp'
    :param key_type: The type of subtask: 'spc', 'pes', or `None`
    :return: A sequence of (task name, task line) pairs
    """
    if task_type not in run_dct:
        return []

    block = run_dct.get(task_type)
    lines = [line.strip() for line in block.strip().splitlines() if line]
    if task_type == "els" or task_type == "proc":
        types = ("spc", "pes")
        assert key_type in types, f"Subtask type {key_type} not in {types}"
        start_key = "ts" if key_type == "pes" else "spc"
        lines = [
            line
            for line in lines
            if line.startswith(start_key) or line.startswith("all")
        ]

    return lines


def task_method_and_basis(
    file_dct: dict[str, str],
    task_key: str,
    key_type: str,
    lvl: str | None = None,
    root: bool = True
) -> tuple[str, str]:
    """Extract the method and basis of a task from a file dictionary

    :param run_dct: The dictionary of a parsed run.dat file
    :param task_key: Task key
    :param subtask_key: Subtask key
    :param key_type: The type of subtask: 'spc', 'pes', or `None`
    :param root: Whether to return the root method/basis (inp level for SP energies)
    :return: A sequence of (task name, task line) pairs
    """
    run_dct = parse_run_dat(file_dct.get("run.dat"))
    thys_dct = parse_theory_dat(file_dct.get("theory.dat"))

    key = "inplvl" if root and task_key in INP_DIR_TASKS else "runlvl"

    if lvl is None:
        task_lines = task_lines_from_run_dict(
            run_dct, task_type="els", key_type=key_type
        )
        task_line = next(line for line in task_lines if task_key in line)
        field_dct = parse_task_fields(task_line)
        lvl = field_dct.get(key)

    thy_dct = thys_dct.get(lvl)
    method = thy_dct.get("method")
    basis = thy_dct.get("basis")
    return (method, basis)


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


def subtask_key_type(key: str) -> str:
    """Determine the type ('pes' or 'spc') of a subtask key.

    :param key: Subtask key
    :return: Subtask key type
    """
    return "spc" if is_species_subtask_key(key) else "pes"


def is_species_subtask_key(key: str) -> bool:
    """Determine whether a subtask key represents a species.

    :param key: Subtask key
    :return: `True` if yes, `False` if no
    """
    return key.strip().isdigit()


def format_subtask_key(key: str) -> str:
    """Format a subtask key

    Examples:
        '1'           ->  '01'          # species
        '1: 2'        ->  '01_2'        # channel
        '1: 1-10'     ->  '01_1-10'     # multi-channel
        '1: 1-10,13'  ->  '01_1-10.13'  # multi-channel
        'all'         ->  'all'         # all

    :param key: The key to parse
    :return: The parsed components
    """
    int_fmt_ = "{:02d}".format

    chn_key = ppc.integer + pp.Suppress(":") + ppc.integer + pp.Suppress(pp.StringEnd())
    pes_key = ppc.integer + pp.Suppress(":") + pp.SkipTo(pp.StringEnd())
    spc_key = ppc.integer + pp.Suppress(pp.StringEnd())
    all_key = pp.Literal(ALL_KEY) + pp.Suppress(pp.StringEnd())
    expr = (chn_key | pes_key | spc_key | all_key) + pp.StringEnd()
    res = expr.parseString(key).as_list()

    # Format integer values
    res = [int_fmt_(x) if isinstance(x, int) else x for x in res]
    assert all(isinstance(x, str) for x in res)

    # Replace special characters
    res = [x.replace(",", ".") for x in res]

    # Join with underscores
    return "_".join(res)


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
    delim = pp.WordEnd() ^ pp.Literal(",")
    #delim = pp.Suppress(pp.White() | ",")
    expr = pp.DelimitedList(entry, delim=delim)
    # expr = pp.OneOrMore(entry + pp.Optional(delim))
    idxs = []
    for res in expr.parseString(inp).as_list():
        if isinstance(res, int):
            idxs.append(res)
        else:
            start, stop = res
            idxs.extend(range(start, stop + 1))
    return list(mit.unique_everseen(idxs))


def format_index_series(idxs: Sequence[int]) -> str:
    r"""Parse a sequence of indices from a string separated by commas and newlines,
    with ranges indicated by 'x-y'

    Example:
        Input: (1, 3, 5, 6, 7, 8, 9, 11, 13, 14, 23, 27, 28, 29)
        Output: '1,3,5-9,11,13-14,23,27-29'
    """
    if not idxs:
        return ""

    def _format_contiguous(idxs_: Sequence[int]) -> str:
        if len(idxs_) > 1:
            start_idx, *_, stop_idx = idxs_
            return f"{start_idx}-{stop_idx}"

        assert len(idxs_) == 1
        (idx,) = idxs_
        return f"{idx}"

    idxs = sorted(idxs)
    counter = itertools.count()
    out = ",".join(
        _format_contiguous(list(g))
        for _, g in itertools.groupby(idxs, lambda idx: idx - next(counter))
    )
    return out


def block_expression(keyword: str, key: str = "content") -> pp.ParseExpression:
    """Parse a block from an AutoMech input file

    :param inp: The input file contents
    :param keyword: The block keyword
    :return: _description_
    """
    start = pp.Keyword(keyword)
    end = pp.Keyword("end") + pp.Keyword(keyword)
    return pp.Suppress(... + start) + pp.SkipTo(end)(key) + pp.Suppress(end)
