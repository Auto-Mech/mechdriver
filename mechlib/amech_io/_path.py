"""
  Handle the generation of necessary paths for various things
"""

import os
import random
import autofile
import automol
from mechanalyzer.inf import spc as sinfo
from mechlib.amech_io import printer as ioprinter


# Set paths to MESS jobs
def messrate_path(prefix, pes_formula, sub_pes_idx):
    """ Build a simple mess path using the run prefix
    """
    pes_str = '{}_{}'.format(pes_formula, sub_pes_idx)
    return os.path.join(prefix, 'MESSRATE', pes_str)


def messpf_path(prefix, inchi):
    """ Build a simple mess path using the run prefix
    """
    spc_formula = automol.inchi.formula_string(inchi)
    ich_key = automol.inchi.inchi_key(inchi)
    return os.path.join(prefix, 'MESSPF', spc_formula, ich_key)


def thermo_paths(spc_dct_i, run_prefix, idx):
    """ Set up the path for saving the pf input and output.
        Placed in a MESSPF, NASA dirs high in run filesys.
    """

    # Get the formula and inchi key
    spc_info = sinfo.from_dct(spc_dct_i)
    spc_formula = automol.inchi.formula_string(spc_info[0])
    ich_key = automol.inchi.inchi_key(spc_info[0])

    # PF
    bld_pf_path = job_path(run_prefix, 'PF', locs_idx=idx)
    spc_pf_path = os.path.join(bld_pf_path, spc_formula, ich_key)

    # NASA polynomials
    bld_nasa_path = job_path(run_prefix, 'NASA', locs_idx=idx)
    spc_nasa_path = os.path.join(bld_nasa_path, spc_formula, ich_key)

    ioprinter.debug_message(
        'preparing thermo for:', spc_info[0],
        bld_pf_path, spc_formula, ich_key)

    return spc_pf_path, spc_nasa_path


def ckin_path():
    """ Set path and make directories
    """
    starting_path = os.getcwd()
    path = os.path.join(starting_path, 'ckin')
    if not os.path.exists(path):
        os.makedirs(path)

    return path


# Helper
def job_path(prefix, job, locs_idx=None, print_path=False):
    """ Create the path for some type of job
    """

    # Initialize the build object
    bld_fs = autofile.fs.build(prefix)

    # Determine the index for the locs if not provided
    if locs_idx is not None:
        assert isinstance(locs_idx, int), (
            'locs idx {} is not an integer'.format(locs_idx)
        )
    else:
        locs_idx = random.randint(0, 9999999)
        bld_fs[-1].create([job, locs_idx])
        # bld_fs[-1].create([job, 0])
        # existing_locs = bld_fs[-1].existing()
        # current_idxs = tuple(idx for [name, idx] in existing_locs
        #                      if name == job)
        # locs_idx = max(current_idxs) + 1

    bld_locs = [job, locs_idx]
    bld_path = bld_fs[-1].path(bld_locs)

    if print_path:
        print('Path for {} Job:'.format(job))
        print(bld_path)

    return bld_path
