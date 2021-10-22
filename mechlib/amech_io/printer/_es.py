"""
ES routines prints
"""

from itertools import chain
import automol
from phydat import phycon
from mechlib.amech_io.printer._print import info_message


def energy(ene):
    """ Print an energy value
    """
    if ene is not None:
        info_message(f"Energy [au]: {ene:<16.8f}", newline=1)


def geometry(geom):
    """ a
    """
    geom_str = automol.geom.string(geom)
    info_message(f'Geometry [Angstrom]:\n{geom_str}', newline=1)


def gradient(grad):
    """ a
    """
    grad_str = automol.util.mat.string(grad)
    info_message(f'Gradient [au]:\n{grad_str}', newline=1)


def frequencies(freqs):
    """ Print out the Harmonic frequencies and ZPVE
    """
    if freqs is not None:
        freq_str = automol.util.vec.string(
            freqs, num_per_row=6, val_format='{0:>12.3f}')
        harm_zpe = (sum(freqs) / 2.0) * phycon.WAVEN2KCAL
        info_message(f'\nHarmonic frequencies [cm-1]:\n{freq_str}')
        info_message(f'\nHarmonic ZPVE [kcal mol-1]: {harm_zpe}')


def molecular_properties(dmom, polar):
    """ a
    """
    if dmom is not None and polar is not None:
        dmom_str = automol.util.vec.string(dmom)
        polar_str = automol.util.mat.string(polar)
        info_message(f'Dipole Moment [Debye]:\n{dmom_str}', newline=1)
        info_message(f'Polarizability []: {polar_str}')


def constraint_dictionary(dct):
    """ a
    """
    if dct is not None:
        info_message('Contraint dictionary:', dct)


def existing_path(label, path):
    """ a
    """
    info_message(label, 'found and saved previously in', path)


def initial_geom_path(label, path):
    """ a
    """
    info_message(label, 'using geom from', path)


def bad_conformer(reason):
    """ a
    """
    info_message(
        f'- Geometry is {reason}. Conformer will not be saved.')


def diverged_ts(param, ref_param, cnf_param):
    """ a
    """
    info_message(
        "- Transition State conformer has",
        "diverged from original structure of",
        f"{param} {ref_param:.3f} with angle {cnf_param:.3f}")


def bad_equil_ts(cnf_dist, equi_bnd):
    """ a
    """
    info_message(
        " - Transition State conformer has",
        "converged to an",
        "equilibrium structure with dist",
        f" {cnf_dist:.3f} comp with equil {equi_bnd:.3f}")


def save_conformer(cnf_save_path):
    """ a
    """
    info_message(" - Geometry is unique. Saving...")
    info_message(f" - Save path: {cnf_save_path}")


def save_conformer_energy(sp_save_path):
    """ a
    """
    info_message(f" - Saving energy of unique geometry at {sp_save_path}...")


def save_symmetry(sym_save_path):
    """ a
    """
    info_message(" - Saving structure in a sym directory at path "
                 f"{sym_save_path}")


def already_running(task, path):
    """ a
    """
    info_message(f'{task} already running in {path}')


def save_reference(save_path):
    """ a
    """
    info_message(" - Saving reference geometry...")
    info_message(f" - Save path: {save_path}")


def run_rotors(run_tors_names, const_names):
    """ a
    """
    info_message(
        'Running hindered rotor scans for the following rotors...',
        newline=1)
    for names in run_tors_names:
        info_message(names)
    if const_names is not None:
        if set(list(chain(*run_tors_names))) == set(const_names):
            info_message(
                'User requested all torsions of system will be fixed.',
                newline=1)


def save_irc(save_path):
    """ a
    """
    info_message(" - Saving IRC...")
    info_message(f" - Save path: {save_path}")


def save_geo(save_path):
    """ a
    """
    info_message(" - Saving geoemetry...")
    info_message(f" - Save path: {save_path}")


def save_energy(sp_save_path):
    """ a
    """
    info_message(" - Saving energy...")
    info_message(f" - Save path: {sp_save_path}")


def save_anharmonicity(geo_save_path):
    """ a
    """
    info_message(" - Saving anharmonicities...")
    info_message(f" - Save path: {geo_save_path}")


def save_frequencies(save_path):
    """ a
    """
    info_message(" - Saving frequencies...")
    info_message(f" - Save path: {save_path}")


def save_gradient(save_path):
    """ a
    """
    info_message(' - Gradient found in Hessian job output.')
    info_message(" - Saving gradient...")
    info_message(f" - Save path: {save_path}")
