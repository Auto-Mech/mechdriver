""" Write NST input
"""

import os
from ioformat import build_mako_str
from nst_io._format import format_geometry
from nst_io._format import format_hessian


# OBTAIN THE PATH TO THE DIRECTORY CONTAINING THE TEMPLATES #
SRC_PATH = os.path.dirname(os.path.realpath(__file__))
TEMPLATE_PATH = os.path.join(SRC_PATH, 'templates')


def input_file(job_type, geo, zero_ene,
               ene_grid=100.0, ene_max=50000.0,
               ang_grid=6, ang_max=300,
               h_so=84.0, qe_scale_factor=1.0,
               cuts=False,
               comment=None):
    """ Writes a string for the input file for NST.

        :param geo: guess geometry for a NST TS
        :type geo: automol.geom object
        :param job_type: NST job to run
        :type job_type: str
        :param geo: geometry of transition state
        :type geo: automol.geom object
        :param zero_ene: energy at infinite separation (in Hartrees)
        :type zero_ene: float
        :param ene_gro: Energy grid spacing (cm-1)
        :param ene_max: maximum grid energy (cm-1)
        :param so_coup: Spin-orbit coupling in cm-1
        :param so_scale: scaling factor for state counts in output
        :param cuts: cuts along vibrational modes (only for HESSC and HESSR)
        :param comment: comment line at top of file
        :rtype: str
    """

    # Format the molecule info
    natoms, geo_str, mass_str = format_geometry(geo)

    # Format flags
    cuts_flag = 'T' if cuts else 'F'

    # Format comment line
    comment = comment if comment is not None else 'NST Run'

    # Create a fill value dictionary
    inp_keys = {
        'comment': comment,
        'job_type': job_type,
        'natoms': natoms,
        'geo_str': geo_str,
        'mass_str': mass_str,
        'zero_ene': zero_ene,
        'ene_grid': ene_grid,
        'ene_max': ene_max,
        'ang_grid': ang_grid,
        'ang_max': ang_max,
        'h_so': h_so,
        'qe_scale_factor': qe_scale_factor,
        'cuts_flag': cuts_flag,
    }

    return build_mako_str(
        template_file_name='input.mako',
        template_src_path=TEMPLATE_PATH,
        template_keys=inp_keys)


def cartesian_hessian_file(hess):
    """ Write the NST cartesian file.
    """
    return format_hessian(hess)
