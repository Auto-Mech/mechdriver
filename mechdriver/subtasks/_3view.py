"""Visualize results."""

from collections.abc import Callable
from pathlib import Path

import automol
import pyparsing as pp
from matplotlib import pyplot
from pyparsing import pyparsing_common as ppc

import autofile

from . import fs


def display(
    task_key: str,
    subtask_key: str,
    path: str | Path = ".",
    runlvl: str | None = None,
    animate: bool = True,
) -> None:
    """Display results for a task.

    :param task_key: Task key
    :param subtask_key: Subtask key
    :param path: Path to AutoMech run directory, containing `inp/` folder
    :param runlvl: Specify the run level (needed if there are duplicate tasks)
    """
    if not task_has_display_function(task_key):
        raise NotImplementedError(f"Display for {task_key} is not yet implemented...")

    disp_ = task_display_function(task_key)

    for task_path in fs.task_paths(task_key, subtask_key, path=path, runlvl=runlvl):
        disp_(task_path)


def _display_init_geom(task_path: str) -> None:
    cnf_fs = autofile.fs.conformer(task_path)
    for cnf_loc in cnf_fs[-1].existing():
        geo = cnf_fs[-1].file.geometry.read(cnf_loc)
        automol.geom.display(geo)


def _display_conf_hess(task_path: str) -> None:
    cnf_fs = autofile.fs.conformer(task_path)
    geo = norm_coo = None
    for cnf_loc in cnf_fs[-1].existing():
        path = cnf_fs[-1].path(cnf_loc)
        print("----")
        print(path)
        geo = cnf_fs[-1].file.geometry.read(cnf_loc)
        hess = cnf_fs[-1].file.hessian.read(cnf_loc)
        freqs, norm_coos = automol.geom.vibrational_analysis(geo, hess)
        norm_coo = norm_coos[:, 0]
        print("Lowest frequency mode:")
        print(f"  frequency: {freqs[0]}")
        print()
    if geo is not None:
        automol.geom.display(geo, mode=norm_coo)


def _display_find_ts(task_path: str) -> None:
    _display_conf_hess(task_path)
    zma_fs = autofile.fs.zmatrix(task_path)
    for zma_loc in zma_fs[-1].existing():
        zma_path = zma_fs[-1].path(zma_loc)
        scan_fs = autofile.fs.scan(zma_path)
        for scan_loc in scan_fs[-2].existing():
            traj = scan_fs[-2].file.trajectory.read(scan_loc)
            if traj:
                geos, comments = zip(*traj)
                enes, coords = zip(*map(parse_scan_trajectory_comment, comments))
                # Plot the energy
                pyplot.plot(coords, enes, marker="o")
                pyplot.show()
                # Show structures
                print(f"Coordinate {coords[0]} geometry: ")
                automol.geom.display(geos[0])
                print(f"Coordinate {coords[len(coords) // 2]} geometry: ")
                automol.geom.display(geos[len(coords) // 2])
                print(f"Coordinate {coords[-1]} geometry: ")
                automol.geom.display(geos[-1])


def _display_rpath_scan(task_path: str) -> None:
    """Display results for a "find_ts" task.

    :param task_path: Task path
    """
    for scan_fs in autofile.fs.iterate_managers(
        task_path, ["CONFORMER", "ZMATRIX"], "SCAN"
    ):
        irc_loc = next(loc for loc in scan_fs[1].existing() if loc[-1][0] == "IRC")
        irc_geos, *_ = zip(*scan_fs[1].file.trajectory.read(irc_loc))
        automol.geom.display_trajectory(irc_geos)


TASK_DISPLAY_FUNCTION = {
    "find_ts": _display_find_ts,
    "conf_hess": _display_conf_hess,
    "conf_opt": _display_conf_hess,
    "rpath_scan": _display_rpath_scan,
    "init_geom": _display_init_geom,
}


def task_has_display_function(task_key: str) -> bool:
    """Determine if task has display function.

    :param task_key: Task key
    :return: `True` if it does, otherwise `False`
    """
    return task_key in TASK_DISPLAY_FUNCTION


def task_display_function(task_key: str) -> Callable[[str], None]:
    """Get the appropriate task display function.

    :param task_key: Task key
    :return: Display function
    """
    return TASK_DISPLAY_FUNCTION.get(task_key)


# Helpers
ENERGY_FIELD = pp.Literal("energy:")
ENERGY_VALUE = ppc.number("energy")
COORD_FIELD = pp.SkipTo(":", include=True)
COORD_VALUE = ppc.number("coord")
SCAN_COMMENT_EXPR = ENERGY_FIELD + ENERGY_VALUE + COORD_FIELD + COORD_VALUE


def parse_scan_trajectory_comment(comment: str):
    """Parse a scan trajectory comment."""
    res = SCAN_COMMENT_EXPR.parse_string(comment)
    return res.get("energy"), res.get("coord")
