"""
Executes the automation part of 1DMin
"""

from routines.trans.lj._geom import get_geometry
from routines.trans.lj._util import util
from routines.trans.runner import runner as ljrunner
from lib import filesys


def onedmin(spc_name, bath_name,
            spc_dct, pf_levels,
            thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix):
    """ Run the task
    """

    spc_info = filesys.inf.get_spc_info(spc_dct[spc_name])
    bath_info = filesys.inf.get_bath_info(spc_dct[bath_name])

    thy_info = ''

    _runlj(spc_info, bath_info,
           spc_dct, pf_levels,
           thy_dct, etrans_keyword_dct,
           run_prefix, save_prefix)
    _savelj(spc_info, bath_info,
            spc_dct, pf_levels,
            thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix)


def _runlj(spc_name, bath_name,
           spc_dct, pf_levels,
           thy_dct, etrans_keyword_dct,
           run_prefix, save_prefix):
    """ Run the Lennard-Jones parameters
    """

    # Set the new information
    chg, mult = util.trans_chg_mult(tgt_info, bath_info)

    njobs = etrans_keyword_dct['njobs']

    # Search save file system for LJ params
    SIGMA, EPSILON = py1dmin.ljroutines.read_lj_from_save(
        target_save_prefix, target_info, theory_info)

    # Run 1DMin to calculate the LJ parameters
    if SIGMA is None and EPSILON is None:  # or NSAMP <= nsamp_run:

        # Obtain the geometry for the target and bath
        target_geo = util.get_geometry(
            tgt_info, thy_info,
            tgt_save_prefix,
            conf=conf)
        bath_geo = py1dmin.ljroutines.get_geometry(
            bath_info, thy_info,
            bath_save_prefix,
            conf=conf)

        # Write the params to the run file system
        FS_THEORY_INFO = [THEORY_INFO[1],
                          THEORY_INFO[2],
                          moldr.util.orbital_restriction(
                              target_info, THEORY_INFO)]
        tgt_run_fs = autofile.fs.species(RUN_PREFIX)
        tgt_run_fs[-1].create(tgt_info)
        tgt_run_path = tgt_run_fs[-1].path(tgt_info)
        etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
        etrans_run_path = etrans_run_fs[-1].path(FS_THEORY_INFO)
        etrans_run_fs[-1].create(FS_THEORY_INFO)

        # Run an instancw of 1DMin for each processor
        for idx in range(njobs):

            # Build run directory
            ljrunner.make_job_dirs(etrans_run_path, idx)

            # Write the input strings
            
            # Write the files

        # Submit the job
        print('\n\nRunning each OneDMin job...')
        ljrunner.submit_job(DRIVE_PATH, etrans_run_path, NJOBS)


def _savelj(spc_name, bath_name,
            spc_dct, pf_levels,
            thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix):
    """ Save the Lennard-Jones parameters
    """

    # Read the lj from all the outputs
    print('\n\nAll OneDMin jobs finished.')
    print('\nReading the Lennard-Jones parameters...')
    SIGMA, EPSILONS = py1dmin.ljroutines.lj_parameters(etrans_run_path)
    if SIGMA is None and EPSILON is None:
        print('\nNo Lennard-Jones parameters found.')
        print('\n\nExiting OneDMin...')
        sys.exit()
    else:
        # Grab the geometries and zero-energies
        print('\nReading the Lennard-Jones Potential Well Geometries...')
        GEOMS = py1dmin.ljroutines.lj_well_geometries(
            etrans_run_path)
        print('\nReading the Zero-Energies...')
        ENES = py1dmin.ljroutines.zero_energies(
            etrans_run_path)

    # Write the params to the save file system
    py1dmin.ljroutines.write_lj_to_save(
        SIGMA, EPSILON,
        TARGET_SAVE_PREFIX, target_info, THEORY_INFO)
    print('\n\nExiting OneDMin...')
