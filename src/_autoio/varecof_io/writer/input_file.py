"""
Writes the global keyword section of a MESS input file
"""

import os
from ioformat import build_mako_str
from phydat import phycon
import elstruct.writer
from varecof_io.writer import _format as vrcformat


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def tst(nsamp_max, nsamp_min, flux_err, pes_size,
        faces=(0,), faces_symm=1,
        ener_grid=(), amom_grid=()):
    """ Writes the tst.inp file for VaReCoF
        :param nsamp_max: maximum number of samples
        :type nsamp_max: int
        :param nsamp_min: minimum number of samples
        :type nsamp_min: int
        :param flux_err: allowed error in flux during sampling
        :type flux_err: float
        :param pes_size: number of PESs in the calculation ?
        :type pes_size: int
        :param faces: inde ?
        :type fascs: list(int)
        :param faces_symm: indices noting symmetry of face
        :type faces_symm: list(int)
        :param list ener_grid:
        :type ener_grid: list(float)
        :param amom_grid:
        :type amom_grid: list(float)
        :rtype: str
    """

    # Set the energy and angular momentum grids
    if not ener_grid:
        ener_grid = [0, 10, 1.05, 179]
    else:
        assert len(ener_grid) == 4
    if not amom_grid:
        amom_grid = [0, 1, 1.10, 40]
    else:
        assert len(amom_grid) == 4
    ener_grid = vrcformat.format_grids_string(ener_grid, 'ener', 'Kelvin')
    amom_grid = vrcformat.format_grids_string(amom_grid, 'amom', 'au')

    # Set the faces
    faces = vrcformat.format_faces_string(faces)

    # Create dictionary to fill template
    tst_keys = {
        'ener_grid': ener_grid,
        'amom_grid': amom_grid,
        'nsamp_max': nsamp_max,
        'nsamp_min': nsamp_min,
        'flux_err': flux_err,
        'pes_size': pes_size,
        'faces': faces,
        'faces_symm': faces_symm
    }

    return build_mako_str(
        template_file_name='tst.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=tst_keys)


def divsur(rdists,
           npivot1,
           npivot2,
           xyz_pivot1,
           xyz_pivot2,
           frame1=(0, 0, 0, 0),
           frame2=(0, 0, 0, 0),
           r2dists=(),
           d1dists=(), d2dists=(),
           t1angs=(), t2angs=(),
           p1angs=(), p2angs=(),
           phi_dependence=False,
           **conditions):
    """ Writes the divsur.inp file for VaReCoF
         that contains info on the dividing surfaces.
        :param rdists: List of temperatures (in Angstrom)
        :type rdists: float
        :param npivot1: number of pivot points on fragment 1
        :type npivot1: int
        :param npivot2: number of pivot points on fragment 2
        :type npivot2: int
        :param xyz_pivot1: xyz of fragment 1 where pivot pts centered
        :type xyz_pivot1: list(float)
        :param xyz_pivot2: xyz of fragment 2 where pivot pts centered
        :type xyz_pivot2: list(float)
        :param frame1: atom idxs for orientation of fragment 1 frame
        :type frame1: list(int)
        :param frame2: atom idxs for orientation of fragment 2 frame
        :type frame2: list(int)
        :param r2dists: cycle coordinates
        :type r2dists: list(float)
        :param d1dists: cycle coordinates
        :type d1dists: list(float)
        :param d2dists: cycle coordinates
        :type d2dists: list(float)
        :param t1angs: cycle coordinates
        :type t1angs: list(float)
        :param t2angs: cycle coordinates
        :type t2angs: list(float)
        :param list p1angs: cycle coordinates
        :type p1angs: list(float)
        :param list p2angs: cycle coordinates
        :type p2angs: list(float)
        :param phi_dependence: signals if frames allow phi angles
        :type phi_dependence: bool
        :param conditions: criteria for cycles
        :type conditions: dict
        :rtype: string
    """

    # Format values strings for the coordinates
    # Function returns the empty string if list is empty
    r1_string = vrcformat.format_values_string(
        'r1', rdists, conv_factor=phycon.ANG2BOHR)
    r2_string = vrcformat.format_values_string(
        'r2', r2dists, conv_factor=phycon.ANG2BOHR)
    d1_string = vrcformat.format_values_string(
        'd1', d1dists, conv_factor=phycon.ANG2BOHR)
    d2_string = vrcformat.format_values_string(
        'd2', d2dists, conv_factor=phycon.ANG2BOHR)
    t1_string = vrcformat.format_values_string(
        't1', t1angs, conv_factor=phycon.RAD2DEG)
    t2_string = vrcformat.format_values_string(
        't2', t2angs, conv_factor=phycon.RAD2DEG)
    p1_string = vrcformat.format_values_string(
        'p1', p1angs, conv_factor=phycon.RAD2DEG)
    p2_string = vrcformat.format_values_string(
        'p2', p2angs, conv_factor=phycon.RAD2DEG)

    # Fromat the frames
    frame1 = ' '.join([str(val) for val in frame1])
    frame2 = ' '.join([str(val) for val in frame2])

    # Write the pivot point coordinates
    idx1 = 1
    idx2 = 1 + npivot1
    if p1angs:
        phi_dependence = True
    pivot_xyz_string1 = vrcformat.format_pivot_xyz_string(
        idx1, npivot1, xyz_pivot1, phi_dependence=phi_dependence)
    pivot_xyz_string2 = vrcformat.format_pivot_xyz_string(
        idx2, npivot2, xyz_pivot2, phi_dependence=phi_dependence)

    # Calculate the number of cycles
    ncycles = 1
    if r2dists:
        ncycles += 1
    if d1dists:
        ncycles += 1
    if d2dists:
        ncycles += 1
    if p1angs:
        ncycles += 1
    if p2angs:
        ncycles += 1
    if t1angs:
        ncycles += 1
    if t2angs:
        ncycles += 1

    # Determine the string of distance cycles
    if d1dists and d2dists:
        dist_coords_string = 'r11 = r1-(d1+d2)/2\n'
        dist_coords_string += 'r21 = r1-(d1+d2)/2\n'
        dist_coords_string += 'r12 = r2-(d1+d2)/2\n'
        dist_coords_string += 'r22 = r2-(d1+d2)/2'
    elif d1dists and not d2dists:
        dist_coords_string = 'r11 = r1-d1/2\n'
        if npivot1 == 2 and npivot2 == 1:
            dist_coords_string += 'r21 = r1-d1/2'
        elif npivot1 == 1 and npivot2 == 2:
            dist_coords_string += 'r12 = r1-d1/2'
    else:
        dist_coords_string = 'r11 = r1'

    # Build string for conditions
    conditions_string = ''
    nconditions = len(conditions)
    if 'delta_r' in conditions:
        alpha = str(conditions['delta_r'])
        conditions_string += f'(r2-r1) < 0.01 + {alpha}\n'
        conditions_string += f'(r2-r1) > 0.01 + {alpha}'

    # Create dictionary to fill template
    divsur_keys = {
        'npivot1': npivot1,
        'npivot2': npivot2,
        'pivot_xyz_string1': pivot_xyz_string1,
        'pivot_xyz_string2': pivot_xyz_string2,
        'frame1': frame1,
        'frame2': frame2,
        'dist_coords_string': dist_coords_string,
        'nconditions': nconditions,
        'conditions_string': conditions_string,
        'ncycles': ncycles,
        'r1_string': r1_string,
        'r2_string': r2_string,
        'd1_string': d1_string,
        'd2_string': d2_string,
        't1_string': t1_string,
        't2_string': t2_string,
        'p1_string': p1_string,
        'p2_string': p2_string,
    }

    return build_mako_str(
        template_file_name='divsur.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=divsur_keys)


def elec_struct(lib_path, base_name, npot,
                dummy_name='dummy_corr_', lib_name='libcorrpot.so',
                exe_name='run.sh',
                geo_ptt='GEOMETRY_HERE', ene_ptt='molpro_energy'):
    """ Writes the electronic structure code input file for VaReCoF
        Currently code only runs with Molpro
        :rtype: string
    """

    # Write the correction potential strings
    pot_path = os.path.join(lib_path, lib_name)
    pot_params_str = ''
    for i in range(npot):
        pot_params_str += f'{base_name+"_corr_":<42s}{"1":<8s}\n'
        pot_params_str += f'{"ParameterInteger":<42s}{i+1:<8d}\n'
    pot_params_str += f'{dummy_name:<42s}{"1":<8s}\n'

    # Set the exe path
    exe_path = os.path.join(lib_path, exe_name)

    # Create dictionary to fill template
    els_keys = {
        'exe_path': exe_path,
        'geo_ptt': geo_ptt,
        'ene_ptt': ene_ptt,
        'base_name': base_name,
        'pot_path': pot_path,
        'pot_params_str': pot_params_str
    }

    return build_mako_str(
        template_file_name='els.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=els_keys)


def structure(geo1, geo2):
    """ Writes the structure input file for VaReCoF
        :param list geo1: geometry of fragment 1
        :param list geo2: geometry of fragment 2
        :rtype: string
    """

    # Determine linearity of molecule
    struct_type1 = vrcformat.determine_struct_type(geo1)
    struct_type2 = vrcformat.determine_struct_type(geo2)

    # Format the coordinates of the geoetry
    natoms1, coords1 = vrcformat.format_coords(geo1)
    natoms2, coords2 = vrcformat.format_coords(geo2)

    # Create dictionary to fill template
    struct_keys = {
        'struct_type1': struct_type1,
        'natoms1': natoms1,
        'coords1': coords1,
        'struct_type2': struct_type2,
        'natoms2': natoms2,
        'coords2': coords2,
    }

    return build_mako_str(
        template_file_name='struct.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=struct_keys)


def molpro_template(ts_info, mod_var_scn_thy_info, inf_sep_ene, cas_kwargs):
    """ Writes a template file for a MOLPRO energy calculation that is
        vrcformatized in the VRC-TST sampling procedure.
    """

    method_dct = {
        'caspt2': 'rs2',
        'caspt2c': 'rs2c',
        'mrcisd_q': 'mrci'
    }
    method = method_dct[mod_var_scn_thy_info[1]]

    # Set the lines for methods
    if 'rs2' in method:
        method_lines = (
            "if (iterations.ge.0) then",
            f"  {{{method},shift=0.25}}",
            f"  molpro_energy = energy + {-1.0*inf_sep_ene}",
            "else",
            "  molpro_energy = 10.0",
            "endif\n"
            # "show[1,e25.15],molpro_energy"
        )
    else:
        method_lines = (
            "if (iterations.ge.0) then",
            "  {mrci}",
            f"  molpro_energy = energd + {-1.0*inf_sep_ene}",
            "else",
            "  molpro_energy = 10.0",
            "endif\n"
            # "show[1,e25.15],molpro_energy"
        )

    # Hacky nonsense to get the correct elstruct string

    # Get the guess into a string
    init_guess_lines = cas_kwargs.get('gen_lines')[1]
    init_guess_lines = tuple(line for line in init_guess_lines
                             if 'molpro_energy' not in line)
    init_guess_str = '\n'.join(init_guess_lines)

    # Write the string using just method lines
    # Then just grab a portion of the string
    cas_kwargs.update({'gen_lines': {3: method_lines}})

    inp_str = elstruct.writer.energy(
        geo='GEOMETRY_HERE',
        charge=ts_info[1],
        mult=ts_info[2],
        method='casscf',
        basis=mod_var_scn_thy_info[2],
        prog=mod_var_scn_thy_info[0],
        orb_type=mod_var_scn_thy_info[3],
        **cas_kwargs
        )
    inp_str = '\n'.join(inp_str.splitlines()[7:])

    tml_str = (
        init_guess_str +
        '\nGEOMETRY_HERE\n' +
        inp_str
    )

    return tml_str


def mc_flux():
    """ Writes the mc_flux.inp file.

        :rtype: string
    """
    return (
        'MultiInputFile          tst.inp\n'
        'OutputFile              mc_flux.out\n'
        'Face                    0\n'
        'ElectronicSurface       0')


def convert():
    """ Writes the convert.inp file.

        :rtype: string
    """
    return 'MultiInputFile    tst.inp\n\n'


def machinefile(machine_dct):
    """ Writes a string that details what machines to run VaReCoF on
        as well as the number of cores to use for each host.

        :param machine_dct: {host name: num cores}
        :type machine_dct: dict[str: int]
        :rtype: str
    """

    machine_file_str = ''
    for name, nproc in machine_dct.items():
        machine_file_str += f'{name}:{nproc}\n'

    return machine_file_str.rstrip()
