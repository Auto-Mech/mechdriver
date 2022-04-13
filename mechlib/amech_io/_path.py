""" Library to build specifically formatted directory paths for
    various calculations conducted by MechDriver.
"""

import os
import random
import autofile
import automol
from mechanalyzer.inf import spc as sinfo


# Set paths to MESS jobs
def rate_paths(pes_dct, run_prefix):
    """ Set up the path for saveing the input and output of
        MESSRATE calculations.

        Creates paths for two different versions of mess

        Run different types of directories (1 PES)
            - fml-base: Standard base rate calculations (use idx)
            - fml-wext: Well-Extended base calculations (use 10*idx)
    """

    rate_path_dct = {}
    for pes_inf in pes_dct:

        pes_fml, pes_idx, subpes_idx = pes_inf
        rate_path_dct[pes_inf] = {}

        for mess_version in ('v1', 'v2'):
            _pes_str = f'{pes_fml}_{str(pes_idx+1)}_{str(subpes_idx+1)}'
            id1 = '{mess_version}-base'
            id2 = '{mess_version}-wext'
            rate_path_dct[pes_inf].update({
                f'base-{mess_version}': job_path(
                    run_prefix, 'MESS', 'RATE', _pes_str,
                    locs_id=id1),
                f'wext-{mess_version}': job_path(
                    run_prefix, 'MESS', 'RATE', _pes_str,
                    locs_id=id2)
            })

    return rate_path_dct


def thermo_paths(spc_dct, spc_locs_dct, spc_mods, run_prefix):
    """ Set up the path for saving the pf input and output.
        Placed in a MESSPF, NASA dirs high in run filesys.
    """

    thm_path_dct = {}
    for spc_name in spc_locs_dct:
        spc_thm_path_dct = {}
        spc_info = sinfo.from_dct(spc_dct[spc_name])
        spc_formula = automol.chi.formula_string(spc_info[0])
        thm_prefix = [spc_formula, automol.chi.inchi_key(spc_info[0])]
        spc_locs_lst = spc_locs_dct[spc_name]
        for sidx, spc_locs in enumerate(spc_locs_lst, start=1):
            spc_mod_thm_path_dct = {}
            for midx, mod in enumerate(spc_mods):
                idx = sidx * 10 + midx
                spc_mod_thm_path_dct[mod] = (
                    job_path(
                        run_prefix, 'MESS', 'PF',
                        thm_prefix, locs_id=idx),
                    job_path(
                        run_prefix, 'THERM', 'NASA',
                        thm_prefix, locs_id=idx)
                )
            spc_mod_thm_path_dct['mod_total'] = (
                job_path(
                    run_prefix, 'MESS', 'PF',
                    thm_prefix, locs_id=sidx),
                job_path(
                    run_prefix, 'THERM', 'NASA',
                    thm_prefix, locs_id=sidx)
            )
            spc_thm_path_dct[tuple(spc_locs)] = spc_mod_thm_path_dct
        spc_thm_path_dct['spc_total'] = (
            job_path(
                run_prefix, 'MESS', 'PF',
                thm_prefix, locs_id=0),
            job_path(
                run_prefix, 'THERM', 'NASA',
                thm_prefix, locs_id=0)
        )
        thm_path_dct[spc_name] = spc_thm_path_dct
    return thm_path_dct


def output_path(dat, make_path=True, print_path=False, prefix=None):
    """ Create the path for sub-directories locatted in the run directory
        where the MechDriver calculation was launched. These sub-directories
        are used to store various useful output from the MechDriver process.

        :param make_path: physically create directory for path during function
        :type make_path: bool
        :param print_path: print the created path to the screen
        :type print_path: bool
        :param prefix: prefix for directory to be built
        :type prefix: str
        :rtype: str
    """

    # Initialize the path
    starting_path = prefix if prefix is not None else os.getcwd()
    path = os.path.join(starting_path, dat)

    # Make and print the path, if requested
    if make_path:
        if not os.path.exists(path):
            os.makedirs(path)
    if print_path:
        print(f'output path for {dat}: {path}')

    return path


def job_path(prefix, prog, job, fml,
             locs_id=None, make_path=True, print_path=False):
    """ Create the path for various types of calculations for
        a given species or PES.

        :param prefix: root prefix to run/save filesyste,
        :type prefix: str
        :param prog: name of the program(s) called in the job
        :type prog: str
        :param fml: stoichiometry of the species/PES associate with job
        :fml type: str
        :param locs_id: some identifier for leaf of the build filesys
        :type locs_id: int
        :param make_path: physically create directory for path during function
        :type make_path: bool
        :param print_path: print the created path to the screen
        :type print_path: bool
        :rtype: str
    """

    # Initialize the build object
    prog_prefix = os.path.join(prefix, prog)
    bld_fs = autofile.fs.build(prog_prefix)

    # Determine the index for the locs if not provided
    if locs_id is None:
        locs_id = str(random.randint(0, 9999999))
    locs_id = str(locs_id)

    if not isinstance(fml, str):
        fml = '-'.join(fml)

    # Build the path
    bld_locs = [job, fml, locs_id]
    bld_path = bld_fs[-1].path(bld_locs)

    # Make and print the path, if requested
    if make_path:
        bld_fs[-1].create([job, fml, locs_id])
    if print_path:
        print(f'Path for {prog}/{job}/{locs_id} Job:')
        print(bld_path)

    return bld_path
