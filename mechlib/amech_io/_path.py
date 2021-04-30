"""
  Handle the generation of necessary paths for various things
"""

import os
import random
import autofile
import automol
from mechanalyzer.inf import spc as sinfo
from mechlib.amech_io import printer as ioprinter


# specialized stuff (to delete)
def write_cwd_pf_file(mess_str, inchi, fname='pf.inp'):
    """ Write a copy of the MESS file in the current working directory
    """

    # Set starting path
    starting_path = os.getcwd()

    # Set the MESS paths and build dirs if needed
    jobdir_mess_path = messpf_path(starting_path, inchi)
    if not os.path.exists(jobdir_mess_path):
        os.makedirs(jobdir_mess_path)

    # Write the files
    file_path = os.path.join(jobdir_mess_path, fname)
    if os.path.exists(file_path):
        for i in range(1, 51):
            if not os.path.exists(file_path+str(i+1)):
                fin_path = file_path+str(i+1)
                break
    else:
        fin_path = file_path
    with open(fin_path, 'w') as file_obj:
        file_obj.write(mess_str)

    ioprinter.saving('MESS input copy', fin_path)

    return fin_path


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
            thm_prefix = [spc_formula, automol.inchi.inchi_key(spc_info[0])]
            print('thm_prefix test:', thm_prefix)
            thm_path[mod] = (
                job_path(run_prefix, 'MESS', 'PF', thm_prefix, locs_idx=idx),
                job_path(run_prefix, 'THERM', 'NASA', thm_prefix, locs_idx=idx)
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
        print('ckin path:'.format(path))
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

    if not isinstance(fml, str):
        fml = '-'.join(fml)
    # Build the path
    print('job path test:', job, fml, locs_idx)
    bld_locs = [job, fml, locs_idx]
    print('job path test 2', bld_locs)
    bld_path = bld_fs[-1].path(bld_locs)
    print('job path test 3', bld_path)

    # Make and print the path if desired
    if make_path:
        bld_fs[-1].create([job, fml, locs_idx])
    if print_path:
        print('Path for {}/{} Job:'.format(prog, job))
        print(bld_path)

    return bld_path
