"""
 functions for reading and writing the filesystem
"""

import automol
import elstruct
import autofile


JOB = elstruct.Job.OPTIMIZATION
RUN_FS = autofile.fs.run('.')
PREFIX = '.'
PROG = 'gaussian09'
METHOD = 'wb97xd'
BASIS = '6-31g*'
ORB_REST = 'R'
THY_LOCS = [METHOD, BASIS, ORB_REST]
CNF_LOCS = ['czW8UnlVcy-aj']
ZMA_LOCS = [0]

OBJECTS = (
    autofile.fs.FileAttributeName.ENERGY,
    autofile.fs.FileAttributeName.GEOM,
)


# CENTRAL SAVING FUNCTION
def save_hr_data(opt_ret, scn_fs, scn_locs):
    """ save info for the hindered rotor
    """

    _save_geom(opt_ret, scn_fs, scn_locs)
    _save_zmatrix(opt_ret, scn_fs, scn_locs)
    _save_energy(opt_ret, scn_fs, scn_locs)


def save_instability(ret, thy_locs, cnf_locs=None, zma_locs=None):
    """ Save instability structure
    """

    # Build filesystem locs and objects
    cnf_locs, zma_locs = _generate_locs(cnf_locs=cnf_locs, zma_locs=zma_locs)

    cnf_fs = autofile.fs.conformer(PREFIX)
    # make paths and pass?
    zma_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))
    sp_fs = autofile.fs.single_point(cnf_fs[-1].path(cnf_locs))
    instab_fs = autofile.fs.zmatrix(cnf_fs[-1].path(cnf_locs))

    # Save data
    _save_geom(ret, cnf_fs, cnf_locs)
    _save_zmatrix(ret, zma_fs, zma_locs)
    # _save_reaction(zrxn)
    _save_energy(ret, sp_fs, thy_locs)
    _save_instab(ret, instab_fs, thy_locs)
