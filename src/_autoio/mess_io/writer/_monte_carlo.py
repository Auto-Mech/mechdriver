"""
Writes MESS input for a monte carlo partition function calculation
"""

import os
from phydat import phycon
import automol.geom
import automol.util.matrix
from ioformat import build_mako_str
from ioformat import remove_trail_whitespace
from ioformat import indent
from mess_io.writer import _format as messformat


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SECTION_PATH = os.path.join(TEMPLATE_PATH, 'sections')
MONTE_CARLO_PATH = os.path.join(SECTION_PATH, 'monte_carlo')


def monte_carlo_species(geo, sym_factor, elec_levels,
                        flux_mode_str, data_file_name,
                        ref_config_file_name='',
                        ground_ene=None, reference_ene=None,
                        freqs=(), excluded_volume_factor=None,
                        use_cm_shift=False):
    """ Writes a monte carlo species section

        :param geo: geometry of species
        :type geo: list
        :param sym_factor: symmetry factor of species
        :type sym_factor: float
        :param elec_levels: energy and degeneracy of atom's electronic states
        :type elec_levels: list(float)
        :param flux_mode_str: MESS-format `FluxionalMode` sections for rotors
        :type flux_mode_str: str
        :param data_file_name: Name of data file with molecular info
        :type data_file_name: str
        :param ground_ene: energy relative to reference n PES (kcal.mol-1)
        :type ground_ene: float
        :param reference_ene: harmonic ZPVE used for MC PF (kcal.mol-1)
        :type reference_ene: float
        :param freqs: vibrational frequencies (cm-1)
        :type freqs: list(float)
        :param use_cm_chift: signal to include a CM shift
        :type use_cm_shift: bool
        :rtype: str
    """

    # Format the molecule specification section
    atom_list = messformat.molec_spec_format(geo)

    # Build a formatted frequencies and elec levels string
    nlevels, levels = messformat.elec_levels_format(elec_levels)
    levels = indent(levels, 2)
    if freqs:
        nfreqs, freqs = messformat.freqs_format(freqs)
    else:
        nfreqs = 0

    # Check if reference config name is present
    # if nfreqs > 0:
    #     assert ref_config_file_name, (
    #         'Must provide a reference configuration file'
    #         'if no Hessians given'
    #     )

    # Indent various strings string if needed
    flux_mode_str = messformat.indent(flux_mode_str, 4)

    # Format the energies
    ref_ene_str = f'{reference_ene:.2f}' if reference_ene is not None else None
    gnd_ene_str = f'{ground_ene:.2f}' if ground_ene is not None else None

    # Create dictionary to fill template
    monte_carlo_keys = {
        'atom_list': atom_list,
        'sym_factor': sym_factor,
        'flux_mode_str': flux_mode_str,
        'data_file_name': data_file_name,
        'ref_config_file_name': ref_config_file_name,
        'reference_ene': ref_ene_str,
        'ground_ene': gnd_ene_str,
        'nlevels': nlevels,
        'levels': levels,
        'nfreqs': nfreqs,
        'freqs': freqs,
        'excluded_volume_factor': excluded_volume_factor,
        'use_cm_shift': use_cm_shift
    }

    return build_mako_str(
        template_file_name='monte_carlo.mako',
        template_src_path=MONTE_CARLO_PATH,
        template_keys=monte_carlo_keys)


def monte_carlo_data(geos, enes, grads=(), hessians=()):
    """ Writes the string for an auxliary data file required for
        Monte Carlo calculations in MESS that contains the
        geoetries, energies, gradients, and Hessians obtained
        from Monte Carlo sampling of the fluxional modes.

        :param geos: geometries from sampling
        :type geos: list
        :param enes: energies from energies
        :type enes: list(float)
        :param grads: gradients from sampling
        :type grads: list
        :param hessians: Hessians from sampling
        :type hessians: list
        :rtype: str
    """

    if not grads and not hessians:
        assert len(geos) == len(enes)
    elif grads or hessians:
        assert grads and hessians
        assert len(geos) == len(enes) == len(grads) == len(hessians)

    dat_str = ''
    for idx, _ in enumerate(geos):
        dat_str += 'Sampling point' + str(idx+1) + '\n'
        dat_str += 'Energy' + '\n'
        dat_str += f'{enes[idx]:.8f}\n'
        dat_str += 'Geometry' + '\n'
        dat_str += messformat.mc_geometry_format(geos[idx]) + '\n'
        if grads:
            grad_str = automol.util.matrix.string(
                grads[idx], val_format='{0:>16.12f}')
            dat_str += 'Gradient'+'\n'
            dat_str += grad_str + '\n'
        if hessians:
            hess_str = automol.util.matrix.string(
                hessians[idx], val_format='{0:>16.12f}')
            dat_str += 'Hessian'+'\n'
            dat_str += hess_str+'\n'

        dat_str += '\n'

    # Format string as needed
    if not grads and not hessians:
        dat_str = remove_trail_whitespace(dat_str)
    dat_str = '\n' + dat_str

    return dat_str


def fluxional_mode(atom_indices, span=6.28319):
    """ Writes the string that defines the `FluxionalMode` section for a
        single fluxional mode (torsion) of a species for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param atom_idxs: idxs of atoms involved in fluxional mode
        :type atom_indices: list(int)
        :param span: range from 0.0 to value that mode was sampled over (rad.)
        :type span: float
        :rtype: str
    """

    # Format the aotm indices string
    atom_indices = messformat.format_flux_mode_indices(atom_indices)

    # Format the span
    span_str = f'{span*phycon.RAD2DEG:.1f}'

    # Create dictionary to fill template
    flux_mode_keys = {
        'atom_indices': atom_indices,
        'span': span_str,
    }

    return build_mako_str(
        template_file_name='fluxional_mode.mako',
        template_src_path=MONTE_CARLO_PATH,
        template_keys=flux_mode_keys)
