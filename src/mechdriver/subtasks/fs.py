"""Get filesystem paths."""

from pathlib import Path

import automol

import autofile

from . import util


def conformer_filesystem_from_task_path(task_path: str):
    """Build a conformer filesystem from a task path.

    :param task_path: Task path
    :return: Conformer filesystem
    """
    return autofile.fs.conformer(task_path)


def task_method_and_basis(
    task_key: str,
    subtask_key: str,
    path: str | Path = ".",
    lvl: str | None = None,
    run: bool = False,
) -> list[Path]:
    """Determine the subtask (species or reaction) method and basis.

    :param subtask_key: Subtask key
    :param path: Path to AutoMech run directory, containing `inp/` folder
    :param run: Whether to do this for the run filesystem, instead of the save
    :return: The species or reaction path
    """
    sub_path = subtask_path(subtask_key, path=path, run=run)
    assert sub_path.exists(), (
        f"Path for {subtask_key} {path} does not exist: {sub_path}"
    )

    # Determine the method and basis
    file_dct = util.read_input_files(path)
    key_type = util.subtask_key_type(subtask_key)
    return util.task_method_and_basis(
        file_dct, task_key, key_type, lvl=lvl, root=False
    )


def task_paths(
    task_key: str,
    subtask_key: str,
    path: str | Path = ".",
    lvl: str | None = None,
    run: bool = False,
) -> list[Path]:
    """Determine the subtask (species or reaction) filesystem path.

    :param subtask_key: Subtask key
    :param path: Path to AutoMech run directory, containing `inp/` folder
    :param run: Whether to do this for the run filesystem, instead of the save
    :return: The species or reaction path
    """
    sub_path = subtask_path(subtask_key, path=path, run=run)
    assert sub_path.exists(), (
        f"Path for {subtask_key} {path} does not exist: {sub_path}"
    )

    # Determine the method and basis
    file_dct = util.read_input_files(path)
    key_type = util.subtask_key_type(subtask_key)
    method, basis = util.task_method_and_basis(
        file_dct, task_key, key_type, lvl=lvl, root=True
    )

    def _eq(val1: str | None, val2: str | None) -> bool:
        if val1 is None or val2 is None:
            return val1 == val2
        return val1.lower() == val2.lower()

    thy_fs = autofile.fs.theory(sub_path)
    thy_loc = next(
        loc
        for loc in thy_fs[-1].existing()
        if _eq(loc[0], method) and _eq(loc[1], basis)
    )
    thy_path = Path(thy_fs[-1].path(thy_loc))

    if util.is_species_subtask_key(subtask_key):
        return [thy_path]

    ts_fs = autofile.fs.transition_state(thy_path)
    ts_locs = ts_fs[-1].existing()
    return list(map(Path, map(ts_fs[-1].path, ts_locs)))


def subtask_path(subtask_key: str, path: str | Path = ".", run: bool = False) -> Path:
    """Determine the subtask (species or reaction) filesystem path.

    :param subtask_key: Subtask key
    :param path: Path to AutoMech run directory, containing `inp/` folder
    :param run: Whether to do this for the run filesystem, instead of the save
    :return: The species or reaction path
    """
    path = Path(path)
    file_dct = util.read_input_files(path)
    spc_df = util.parse_species_csv(file_dct.get("species.csv"))
    rxn_dct = util.parse_mechanism_dat(file_dct.get("mechanism.dat"))
    run_dct = util.parse_run_dat(file_dct.get("run.dat"))
    save_path, run_path = util.filesystem_paths_from_run_dict(run_dct, cwd=path)
    root_path = run_path if run else save_path

    if util.is_species_subtask_key(subtask_key):
        spc_idx = int(subtask_key) - 1
        spc_locs = list(spc_df[["inchi", "charge", "mult"]].iloc[spc_idx])
        spc_locs[0] = automol.chi.canonical_enantiomer(spc_locs[0])
        spc_fs = autofile.fs.species(root_path)
        return Path(spc_fs[-1].path(spc_locs))

    assert subtask_key in rxn_dct, f"Reaction not found: {subtask_key}"
    chi_dct = dict(zip(spc_df["name"], spc_df["inchi"], strict=True))
    mul_dct = dict(zip(spc_df["name"], spc_df["mult"], strict=True))
    chg_dct = dict(zip(spc_df["name"], spc_df["charge"], strict=True))

    rxn_chis = [list(map(chi_dct.get, rs)) for rs in rxn_dct.get(subtask_key)]
    rxn_chis = automol.chi.canonical_enantiomer_reaction(*rxn_chis)
    rxn_chgs = [list(map(chg_dct.get, rs)) for rs in rxn_dct.get(subtask_key)]
    rxn_muls = [list(map(mul_dct.get, rs)) for rs in rxn_dct.get(subtask_key)]
    rxn_chis, rxn_chgs, rxn_muls = autofile.schema.loc_maps.sort_together(
        rxn_chis, rxn_chgs, rxn_muls
    )
    ts_mul = util.ts_multiplicity(*rxn_muls)
    rxn_locs = [rxn_chis, rxn_chgs, rxn_muls, ts_mul]
    rxn_fs = autofile.fs.reaction(root_path)
    return Path(rxn_fs[-1].path(rxn_locs))
