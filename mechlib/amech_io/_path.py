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
def thermo_paths(spc_dct, spc_queue, run_prefix):
    """ Set up the path for saving the pf input and output.
        Placed in a MESSPF, NASA dirs high in run filesys.
    """

    thm_paths = []
    for spc_name, (_, mods, _, _) in spc_queue:
        thm_path = {}
        for idx, mod in enumerate(mods):
            spc_info = sinfo.from_dct(spc_dct[spc_name])
            spc_formula = automol.inchi.formula_string(spc_info[0])
            thm_path[mod] = (
                job_path(run_prefix, 'MESS', 'PF', spc_formula, locs_idx=idx),
                job_path(run_prefix, 'THERM', 'NASA', spc_fourmula, locs_idx=idx)
            )
        thm_paths.append(thm_path)

    return thm_paths


def output_path(dat, make_path=True, print_path=False):
    """ Set path and make directories for making additional
        files in the directory where mechdriver drops are submitted
        and its output made
    """

    # Initialize the path
    starting_path = os.getcwd()
    path = os.path.join(starting_path, dat)

    # Make and print the path if desired
    if make_path:
        if not os.path.exists(path):
            os.makedirs(path)
    if print_path:
        print('ckin path:'.format(path)
        print(bld_path)

    return path


def job_path(prefix, prog, job, fml,
             locs_idx=None, make_path=True, print_path=False):
    """ Create the path for some type of job
    """

    # Initialize the build object
    prog_prefix = os.path.join(prefix, prog)
    bld_fs = autofile.fs.build(prog_prefix)

    # Determine the index for the locs if not provided
    if locs_idx is not None:
        assert isinstance(locs_idx, int), (
            'locs idx {} is not an integer'.format(locs_idx)
        )
    else:
        locs_idx = random.randint(0, 9999999)

    # Build the path
    bld_locs = [job, fml, locs_idx]
    bld_path = bld_fs[-1].path(bld_locs)
    
    # Make and print the path if desired
    if make_path:
        bld_fs[-1].create([job, fml, locs_idx])
    if print_path:
        print('Path for {}/{} Job:'.format(prog, job))
        print(bld_path)

    return bld_path
