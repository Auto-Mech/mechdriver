"""
Writes MESS input for a molecule
"""

import os
import automol._deprecated
from ioformat import build_mako_str
from ioformat import indent
from mess_io.writer import _format as messformat


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')
SPECIES_PATH = os.path.join(TEMPLATE_PATH, 'species')
SPEC_INFO_PATH = os.path.join(SPECIES_PATH, 'info')


def core_rigidrotor(
        geo, sym_factor,
        xmat=(), freqs=(), rovib_coups=(),
        rot_dists=(), interp_emax=None):
    """ Writes the string that defines the 'Core' section for a
        rigid-rotor model of a species for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param geo: geometry of species
        :type geo: list
        :param sym_factor: symmetry factor of species
        :type sym_factor: float
        :param interp_emax: max energy to calculate num. of states (kcal.mol-1)
        :type interp_emax: float
        :rtype: str
    """

    # Format the geometry section
    natom, geo = messformat.geometry_format(geo)

    if xmat:
        anharm = messformat.format_xmat(xmat)
        interp_emax = 500
        nfreqs, freqs = messformat.freqs_format(freqs)
    else:
        anharm = ''
        freqs = ''
        nfreqs = ''
    if rovib_coups:
        rovib_coups = messformat.format_rovib_coups(rovib_coups)
        interp_emax = 500
    else:
        rovib_coups = ''
    if rot_dists:
        rot_dists = messformat.format_rot_dist_consts(rot_dists)
    else:
        rot_dists = ''
    # Create dictionary to fill template
    core_keys = {
        'sym_factor': sym_factor,
        'natom': natom,
        'geo': geo,
        'interp_emax': interp_emax,
        'anharm': anharm,
        'nfreqs': nfreqs,
        'freqs': freqs,
        'rovib_coups': rovib_coups,
        'rot_dists': rot_dists,
    }

    return build_mako_str(
        template_file_name='core_rigidrotor.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=core_keys)


def core_multirotor(geo, sym_factor, pot_surf_file, int_rot_str,
                    interp_emax=100, quant_lvl_emax=9):
    """ Writes the string that defines the `Core` section for a
        multidimensional rotor model of a species for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param geo: geometry of species
        :type geo: list
        :param sym_factor: symmetry factor of species
        :type sym_factor: float
        :param pot_surf_file: name of file with PES along rotor (kcal.mol-1)
        :type pot_sur_file: str
        :param int_rot_str: MESS-format strings that define internal rotors
        :type int_rot_str: str
        :param interp_emax: max energy to calculate density/number of states
        :type interp_emax: float
        :param quant_lvl_emax: max energy to calculate quantum energy levels
        :type quant_lvl_emax: float
        :rtype: str
    """

    # Format the geometry section
    natom, geo = messformat.geometry_format(geo)

    # Indent the internal rotor string
    int_rot_str = indent(int_rot_str, 2)

    # Create dictionary to fill template
    core_keys = {
        'sym_factor': sym_factor,
        'natom': natom,
        'geo': geo,
        'pot_surf_file': pot_surf_file,
        'int_rot': int_rot_str,
        'interp_emax': interp_emax,
        'quant_lvl_emax': quant_lvl_emax
    }

    return build_mako_str(
        template_file_name='core_multirotor.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=core_keys)


def core_phasespace(geo1, geo2, sym_factor, stoich,
                    pot_prefactor=10.0, pot_exp=6.0, tstlvl='e'):
    """ Writes the string that defines the `Core` section for a
        phase space theory model of a transition state for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param geo1: geometry of the dissociation species 1
        :type geo1: list
        :param geo2: geometry of the dissociation species 2
        :type geo2: list
        :param sym_factor: symmetry factor of transition state
        :type sym_factor: float
        :param stoich: combined stoichiometry of dissociation species 1 and 2
        :type stoich: str
        :param pot_prefator: factor C0 in potential expression V = -C0/R^n (au)
        :type pot_prefactor: float
        :param pot_exp: power n in potential expression V = -C0/R^n (au)
        :type pot_exp: float
        :param tstlvl: level to resolve the rate constant
        :type tstlvl: str
        :rtype: str
    """

    assert tstlvl in ('e', 'ej', 't')

    # Format the geometry section of each fragment
    natom1, geo1 = messformat.geometry_format(geo1)
    natom2, geo2 = messformat.geometry_format(geo2)

    # Indent the geometry strings
    geo1 = indent(geo1, 2)
    geo2 = indent(geo2, 2)

    # Format the tstlvl string
    assert tstlvl in ('e', 'ej', 't')
    tstlvl = tstlvl.upper()

    # Create dictionary to fill template
    core_keys = {
        'sym_factor': sym_factor,
        'natom1': natom1,
        'geo1': geo1,
        'natom2': natom2,
        'geo2': geo2,
        'stoich': stoich,
        'pot_prefactor': pot_prefactor,
        'pot_exp': pot_exp,
        'tstlvl': tstlvl
    }

    return build_mako_str(
        template_file_name='core_phasespace.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=core_keys)


def core_rotd(sym_factor, flux_file_name, stoich):
    """ Writes the string that defines the `Core` section for a
        variational reaction-coordinate transition-state theory model of a
        transition state for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param sym_factor: symmetry factor of transition state
        :type sym_factor: float
        :param flux_file_name:
        :type flux_file_name: str
        :param stoich: combined stoichiometry of dissociation species 1 and 2
        :type stoich: str
        :rtype: str
    """

    # Create dictionary to fill template
    core_keys = {
        'sym_factor': sym_factor,
        'flux_file_name': flux_file_name,
        'stoich': stoich
    }

    return build_mako_str(
        template_file_name='core_rotd.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=core_keys)


def rotor_hindered(group, axis, symmetry, potential,
                   hmin=None, hmax=None,
                   lvl_ene_max=None,
                   therm_pow_max=None,
                   geo=None,
                   rotor_id='',
                   potential_form='spline'):
    """ Writes the string that defines the `Rotor` section for a
        single hindered rotor of a species for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param group: idxs for the atoms of one of the rotational groups
        :type group: list(int)
        :param axis: idxs for the atoms that make up the rotational axis
        :type axis: list(int)
        :param symmetry: overall symmetry of the torsional motion (potential)
        :type symmetry: int
        :param hmin: minimum value for quantum phase space dimension
        :type hmin: int
        :param hmax: maximum value for quantum phase space dimension
        :type hmax: int
        :param potential: value of the potential along torsion (kcal.mol-1)
        :type potential: list(float)
        :param therm_pow_max: max exp't power in Boltzmann weight
        :type therm_pow_max: int
        :param geo: geometry of the species the rotor exists for
        :type geo: list
        :param rotor_id: name associated with the rotor
        :type rotor_id: str
        :param potential_form: expression the potential should be fit to
        :type potential_form: str
        :rtype: str
    """

    # Format the rotor sections
    fmtd_group = messformat.format_rotor_key_defs(group)
    fmtd_axis = messformat.format_rotor_key_defs(axis)
    npot, fmtd_coords, fmtd_enes = messformat.format_rotor_potential(
        potential)

    # Format the geom
    natom = 1
    if geo is not None:
        natom, geo = messformat.geometry_format(geo)
        geo = indent(geo, 4)

    # Create dictionary to fill template
    rotor_keys = {
        'group': fmtd_group,
        'axis': fmtd_axis,
        'symmetry': symmetry,
        'npotential': npot,
        'pot_coords': fmtd_coords,
        'pot_enes': fmtd_enes,
        'potential_form': potential_form,
        'hmin': hmin,
        'hmax': hmax,
        'lvl_ene_max': lvl_ene_max,
        'therm_pow_max': therm_pow_max,
        'natom': natom,
        'geo': geo,
        'rotor_id': rotor_id
    }

    return build_mako_str(
        template_file_name='rotor_hindered.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=rotor_keys)


def rotor_internal(group, axis, symmetry, grid_size, mass_exp_size,
                   pot_exp_size=5, hmin=13, hmax=101,
                   geo=None, rotor_id=''):
    """ Writes the string that defines the `Rotor` section for a
        single internal rotor of a species for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param group: idxs for the atoms of one of the rotational groups
        :type group: list(int)
        :param axis: idxs for the atoms that make up the rotational axis
        :type axis: list(int)
        :param symmetry: overall symmetry of the torsional motion (potential)
        :type symmetry: int
        :param grid_size: grid_size for statistical weight calculation
        :type grid_size: int
        :param mass_exp_size: num. mass expansion Fourier harmonics
        :type mass_exp_size: int
        :param pot_exp_size: num. potential expansion Fourier harmonics
        :type pot_exp_size: int
        :param hmin: minimum value for quantum phase space dimension
        :type hmin: int
        :param hmax: maximum value for quantum phase space dimension
        :type hmax: int
        :param geo: geometry of the species the rotor exists for
        :type geo: list
        :param rotor_id: name associated with the rotor
        :type rotor_id: str
        :rtype: str
    """

    assert mass_exp_size > 0 and mass_exp_size % 2 == 1, (
        f'Mass exponent size: {mass_exp_size} is not an odd number'
    )
    assert pot_exp_size > 0 and pot_exp_size % 2 == 1, (
        f'Potential exponent size: {pot_exp_size} is not an odd number'
    )

    # Format the sections
    rotor_group = messformat.format_rotor_key_defs(group)
    rotor_axis = messformat.format_rotor_key_defs(axis)

    # Format the geom
    if geo is not None:
        natom, geo = messformat.geometry_format(geo)
        geo = indent(geo, 4)
    else:
        natom = None

    # Create dictionary to fill template
    rotor_keys = {
        'group': rotor_group,
        'axis': rotor_axis,
        'symmetry': symmetry,
        'mass_exp_size': mass_exp_size,
        'pot_exp_size': pot_exp_size,
        'hmin': hmin,
        'hmax': hmax,
        'grid_size': grid_size,
        'natom': natom,
        'geo': geo,
        'rotor_id': rotor_id
    }

    return build_mako_str(
        template_file_name='rotor_internal.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=rotor_keys)


def mdhr_data(pot_dct, freqs=None, nrot=0):
    """ Writes the string for an auxiliary data file for MESS containing
        potentials and vibrational frequencies of a
        multidimensional hindered rotor, up to four dimensions.

        :param pots: potential values along torsional modes of rotor
        :type pots: list(list(float))
        :param freqs: vibrational frequenciess along torsional modes of rotor
        :type freqs: list(list(float))
        :rtype: str
    """

    assert pot_dct, 'Potential has no values'
    # Remap potential so that keys are indices, not vcoord valyes
    pot_obj = automol.data.potent.from_dict(pot_dct)
    pots_byidx = automol.data.potent.dict_(pot_obj, index=True, drop_null=True)
    pots_bykey = automol.data.potent.dict_(pot_obj, index=False, drop_null=True)
    pot_idxs = tuple(pots_byidx.keys())
    key_idx_dct = {
        idx: key for key, idx in zip(pots_bykey.keys(), pots_byidx.keys())
        if pots_byidx[idx] == pots_bykey[key]}
    # Get the dimensions of the MDHR
    ndims = len(pot_idxs[0])
    assert ndims in (1, 2, 3, 4), 'Rotor must have dimension 1-4'

    # Get the number of terms in each rotor of MDHR
    # Basically finds number of terms for each position of potential grid,
    # for (m, n, ...)->dims=(unique vals in pos m, uniquevals in pos n, ...)
    dims = tuple()
    for dim in range(ndims):
        dims += (max((x[dim] for x in pot_idxs))+1,)

    # Get the number of freqs
    if freqs is not None:
        nfreqs = len(list(freqs.values())[0])
        nfreqs -= nrot
    else:
        nfreqs = 0

    # Write top line string with number of points in potential
    pot_head = 'V(kcal/mol)'
    lbls = ['i', 'j', 'k', 'l']
    if ndims == 1:
        num_str = f'{dims[0]:>6d}\n'
        head_str = f'{lbls[0]:>6s}{pot_head:>15s}'
    elif ndims == 2:
        num_str = f'{dims[0]:>6d}{dims[1]:>6d}\n'
        head_str = f'{lbls[0]:>6s}{lbls[1]:>6s}{pot_head:>15s}'
    elif ndims == 3:
        num_str = f'{dims[0]:>6d}{dims[1]:>6d}{dims[2]:>6d}\n'.format(*dims)
        head_str = f'{lbls[0]:>6s}{lbls[1]:>6s}{lbls[2]:>6s}{pot_head:>15s}'
    elif ndims == 4:
        num_str = f'{dims[0]:>6d}{dims[1]:>6d}{dims[2]:>6d}{dims[3]:>6d}\n'
        head_str = (f'{lbls[0]:>6s}{lbls[1]:>6s}'
                    f'{lbls[2]:>6s}{lbls[3]:>6s}{pot_head:>15s}')

    # Add the nofreq line
    if nfreqs > 0:
        freq_str = ' '.join(f'{i+1:d}' for i in range(nfreqs)) + '\n'
        head_str += 'Freqs(cm-1)' + '\n'
    else:
        freq_str = '\n'
        head_str += '\n'

    # Build the lines for each point on the potential
    dat_str = num_str + freq_str + head_str
    for idxs, val in pots_byidx.items():

        if val is not None:

            # Add the idxs for the rotors
            for idx in idxs:
                dat_str += f'{idx+1:>6d}'

            # Add the potential value
            dat_str += f'{val:>15f}'

            # Add any frequencies if necessary
            if freqs is not None:
                if key_idx_dct[idxs] in freqs:
                    for freq in freqs[key_idx_dct[idxs]]:
                        dat_str += f'{freq:>8.1f}'

            dat_str += '\n'

    return dat_str


def umbrella_mode(group, plane, ref_atom, potential,
                  geo=None):
    """ Writes the string that defines the `Umbrella` section for a
        single umbrella mode of a species for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param group: idxs for the atoms of ?
        :type group: list(int)
        :param axis: idxs for the atoms that ?
        :type axis: list(int)
        :param geo: geometry of the species the umbrella mode exists for
        :type geo: list
        :rtype: str
    """

    # Format the sections
    umbr_group = messformat.format_rotor_key_defs(group)
    umbr_plane = messformat.format_rotor_key_defs(plane)
    umbr_npotential, _, umbr_potential = messformat.format_rotor_potential(
        potential)
    ref_atom += 1

    # Format the geom
    if geo is not None:
        natom, geo = messformat.geometry_format(geo)
        geo = indent(geo, 4)
    else:
        natom = None

    # Create dictionary to fill template
    umbr_keys = {
        'group': umbr_group,
        'axis': umbr_plane,
        'ref_atom': ref_atom,
        'npotential': umbr_npotential,
        'potential': umbr_potential,
        'natom': natom,
        'geo': geo,
    }

    return build_mako_str(
        template_file_name='umbrella_mode.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=umbr_keys)


def tunnel_eckart(imag_freq, well_depth1, well_depth2):
    """ Writes the string that defines the 'Tunneling' section for a
        Eckart tunneling model for a transition state for a MESS input file by
        formatting input information into strings a filling Mako template.

        :param imag_freq: imaginary frequency of the TS
        :type imag_freq: float
        :param well_depth1: energy difference: E[TS] - E[reactant well]
        :type well_depth1: float
        :param well_depth2: energy difference: E[TS] - E[product well]
        :type well_depth2: float
        :rtype: str
    """

    # Format the imaginary frequency and well-depth values
    imag_freq = f'{imag_freq:<8.0f}'
    well_depth1 = f'{well_depth1:<8.2f}'
    well_depth2 = f'{well_depth2:<8.2f}'

    # Create dictionary to fill template
    tunnel_keys = {
        'imag_freq': imag_freq,
        'well_depth1': well_depth1,
        'well_depth2': well_depth2
    }

    return build_mako_str(
        template_file_name='tunnel_eckart.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=tunnel_keys)


def tunnel_read(imag_freq, tunnel_file, cutoff_energy=2500.0):
    """ Writes the string that defines the 'Tunneling' section for a
        small curvature tunneling model for a transition state
        for a MESS input file by formatting input information into
        strings a filling Mako template. Currently requires an additional
        auxiliary data file for transmission probabilities generated by
        ProjRot code.

        :param imag_freq: imaginary frequency of the TS
        :type imag_freq: float
        :param tunnel_file: name of data file with transmission probabilities
        :type tunnel_file: str
        :param cutoff_energy: energy to include tunneling density (kcal.mol-1)
        :type cutoff_energy: float
        :rtype: str
    """

    # Format the imaginary frequency value
    imag_freq = f'{imag_freq:<8.0f}'

    # Create dictionary to fill template
    tunnel_keys = {
        'imag_freq': imag_freq,
        'cutoff_energy': cutoff_energy,
        'tunnel_file': tunnel_file
    }

    return build_mako_str(
        template_file_name='tunnel_read.mako',
        template_src_path=SPEC_INFO_PATH,
        template_keys=tunnel_keys)
