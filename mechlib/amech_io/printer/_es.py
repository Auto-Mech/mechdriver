"""
ES routines prints
"""

from itertools import chain
# import logging

import automol
from phydat import phycon
from mechlib.amech_io.printer._tsk import keyword_list
from mechlib.amech_io.printer import info_message

def es_task_header(tsk, spc_name, spc_smiles, es_keyword_dct, thy_dct):
    """ Print the top of the section for an electronic structure task
    """
    print()
    # printer.obj('line_dash', prechar='\n')
    # task_str = (
    #     'Task'
    #     '\n'
    #     '  {}'.format(tsk)
    #     '\n\n'
    #     'Species'
    #     '  Name: {}'.format(spc_name)
    #     '  SMILES: {}'.format(spc_smiles)
    #     '\n\n'
    #     'Options'
    #     keyword_list(keyword_dct, thy_dct=thy_dct)  # need indent 
    # )
    print(task_str)


# Results
def energy(ene):
    """ Print an energy value
    """
    if ene is not None:
        info_message("Energy [au]: {0:<16.8f}".format(ene), newline=1)


def geometry(geom):
    """ a
    """
    geom_str = automol.geom.string(geom)
    info_message('Geometry [Angstrom]:\n{}'.format(geom_str), newline=1)


def gradient(grad):
    """ a
    """
    grad_str = automol.cart.mat.string(grad)
    info_message('Gradient [au]:\n{}'.format(grad_str), newline=1)


def frequencies(freqs):
    """ Print out the Harmonic frequencies and ZPVE
    """
    if freqs is not None:
        freq_str = automol.cart.vec.string(freqs, num_per_row=6)
        harm_zpe = (sum(freqs) / 2.0) * phycon.WAVEN2EH
        info_message('Harmonic frequencies [cm-1]:\n{}'.format(freq_str), newline=1)
        info_message('Harmonic ZPVE [cm-1]: {}'.format(harm_zpe), newline=1)


def molecular_properties(dmom, polar):
    """ a
    """
    if dmom is not None and polar is not None:
        dmom_str = automol.cart.vec.string(dmom)
        polar_str = automol.cart.mat.string(polar)
        info_message('Dipole Moment [Debye]:\n{}'.format(dmom_str), newline=1)
        info_message('Polarizability []: {}'.format(polar_str))


def constraint_dictionary(dct):
    """ a
    """
    if dct is not None:
        info_message('Contraint dictionary: ', dct)


def existing_path(label, path):
    """ a
    """
    info_message(label + ' found and saved previously in ', path)


def initial_geom_path(label, path):
    """ a
    """
    info_message(label + ' using geom from ', path)


def bad_conformer(reason):
    """ a
    """
    info_message(
        '- Geometry is {}. Conformer will not be saved.'.format(reason))


def diverged_ts(param, ref_param, cnf_param):
    """ a
    """
    info_message(
        "- Transition State conformer has",
        "diverged from original structure of",
        "{} {:.3f} with angle {:.3f}".format(
            param, ref_param, cnf_param))


def bad_equil_ts(cnf_dist, equi_bnd):
    """ a
    """
    info_message(
        " - Transition State conformer has",
        "converged to an",
        "equilibrium structure with dist",
        " {:.3f} comp with equil {:.3f}".format(
            cnf_dist, equi_bnd))


def save_conformer(cnf_save_path):
    """ a
    """
    info_message(" - Geometry is unique. Saving...")
    info_message(" - Save path: {}".format(cnf_save_path))


def save_conformer_energy(sp_save_path):
    """ a
    """
    info_message(" - Saving energy of unique geometry at {}...", sp_save_path)


def save_symmetry(sym_save_path):
    """ a
    """
    info_message(" - Saving structure in a sym directory at path {}".format(
        sym_save_path))


def already_running(task, path):
    """ a
    """
    info_message('{} already running in {}'.format(task, path))


def save_reference(save_path):
    """ a
    """
    info_message(" - Saving reference geometryg...")
    info_message(" - Save path: {}".format(save_path))


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
    info_message(" - Save path: {}".format(save_path))


def save_geo(save_path):
    """ a
    """
    info_message(" - Saving geoemetry...")
    info_message(" - Save path: {}".format(save_path))


def save_energy(sp_save_path):
    """ a
    """
    info_message(" - Saving energy...")
    info_message(" - Save path: {}".format(sp_save_path))


def save_anharmonicity(geo_save_path):
    """ a
    """
    info_message(" - Saving anharmonicities...")
    info_message(" - Save path: {}".format(geo_save_path))


def save_frequencies(save_path):
    """ a
    """
    info_message(" - Saving frequencies...")
    info_message(" - Save path: {}".format(save_path))


def save_gradient(save_path):
    """ a
    """
    info_message(' - Gradient found in Hessian job output.')
    info_message(" - Saving gradient...")
    info_message(" - Save path: {}".format(save_path))
