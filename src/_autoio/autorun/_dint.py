""" autorun functions to do the dintmom.pl script

    STEPS:
        (1) optimize geometry
        (2) sample initial coords of the unimol spc
        (3) lj params
        (4) run collid traj
        (5) calc moms

    Auxiliary input:
        from PIPPy: coef.dat, basis.dat
        made by hand?: tinker files
"""

# import dint_io
# import onedmin_io
from autorun._run import from_input_string


INPUT_NAME = 'input'
OUTPUT_NAMES = ('output',)
OPT_OUTPUT_NAMES = ('dint.geo',)


# autorun functions
# def optimized_geom(script_str, run_dir,
#                    geo, basis_str, coef_str):
#     """ Optimize the geometry using DiNT.
#     """
#
#     # input_str = dint_io.writer.optimization_input(geo)
#     input_str = ''
#     aux_dct = {
#         'basis.dat': basis_str,
#         'coef.dat': coef_str
#         # 'tinker.xyz': tinker_xyz
#     }
#
#     output_strs = direct(
#         script_str, run_dir, input_str,
#         aux_dct=aux_dct,
#         input_name='input',
#         output_names=('dint.geo',))
#
#     # Parse out the info
#     opt_geo = dint_io.reader.geo(output_strs[0])
#     rot_consts = dint_io.reader.geo(output_strs[0])
#     energy = dint_io.reader.geo(output_strs[0])
#
#     return opt_geo, rot_consts, energy


def sampling():
    """ Sample the initial coordinates.
        run optimized geometries for severla geoms
        find a dint.brot file or build from fort.80 with brot.x
    """
    return NotImplementedError


def collistion_trajectory():
    """ Trajectory simulations
    """
    return NotImplementedError


def moments():
    """ Moments
        # Write param.inc
        # compile mom.x file to get moments
        maybe combine with collsion trajectory calculation
    """
    return NotImplementedError


def direct(script_str, run_dir, input_str, aux_dct=None,
           input_name=INPUT_NAME,
           output_names=OUTPUT_NAMES):
    """
        :param input_writer: elstruct writer module function for desired job
        :type input_writer: elstruct function
        :param aux_dct: auxiliary input strings dict[name: string]
        :type aux_dct: dict[str: str]
        :param script_str: string of bash script that contains
            execution instructions electronic structure job
        :type script_str: str
        :param run_dir: name of directory to run electronic structure job
        :type run_dir: str
    """
    return from_input_string(
        script_str, run_dir, input_str,
        aux_dct=aux_dct,
        input_name=input_name,
        output_names=output_names)
