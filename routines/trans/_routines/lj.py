"""
Executes the automation part of 1DMin
"""

from routines.trans.lj._geom import get_geometry
from routines.trans.lj._util import util
from routines.trans.runner import runner as ljrunner
from lib import filesys


def onedmin(spc_name, bath_name,
            spc_dct, thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix):
    """ Run the task
    """

    tgt_dct, bath_dct = spc_dct[spc_name], spc_dct[bath_name]
    tgt_info = filesys.inf.get_spc_info(tgt_dct)
    bath_info = filesys.inf.get_spc_info(bath_dct)
    lj_info = filesys.inf.combine_spc_info(tgt_info, bath_info)

    # Build the modified thy objs
    thy_info = filesys.inf.get_es_info(
        etrans_keyword_dct['runlvl'], thy_dct)
    mod_tgt_thy_info = filesys.inf.modify_orb_restrict(tgt_info, thy_info)
    mod_lj_thy_info = filesys.inf.modify_orb_restrict(lj_info, thy_info)

    # Build the conformer filesystem objects
    _, thy_run_path = filesys.build.spc_thy_fs_from_root(
        run_prefix, tgt_info, mod_tgt_thy_info)
    _, thy_save_path = filesys.build.spc_thy_fs_from_root(
        save_prefix, tgt_info, mod_tgt_thy_info)
    cnf_run_fs, _ = filesys.build.cnf_fs_from_prefix(
        thy_run_path, mod_tgt_thy_info, cnf=None)
    cnf_save_fs, cnf_save_locs = filesys.build.cnf_fs_from_prefix(
        thy_save_path, mod_tgt_thy_info, cnf='min')
    cnf_run_paths = filesys.build.cnf_paths_from_locs(
        cnf_run_fs, cnf_save_locs)
    cnf_save_paths = filesys.build.cnf_paths_from_locs(
        cnf_save_fs, cnf_save_locs)

    # Build the energy transfer filesystem objects
    etrans_run_fs, _ = filesys.build.etrans_fs_from_prefix(
        cnf_run_paths[0], bath_info, mod_lj_thy_info)
    etrans_save_fs, etrans_locs = filesys.build.etrans_fs_from_prefix(
        cnf_save_paths[0], bath_info, mod_lj_thy_info)

    # Calculate and save the Lennard-Jones parameters
    _runlj(etrans_run_fs, etrans_save_fs, etrans_locs,
           etrans_keyword_dct)
    # _savelj()


def _runlj(etrans_run_fs, etrans_save_fs, etrans_locs,
           etrans_keyword_dct):
    """ Run the Lennard-Jones parameters
    """

    # Pull stuff from dct
    njobs = etrans_keyword_dct['njobs']
    overwrite = etrans_keyword_dct['overwrite']

    eps_exists = etrans_save_fs[-1].file.epsilon.exists(etrans_locs)
    sig_exists = etrans_save_fs[-1].file.sigma.exists(etrans_locs)
    if not eps_exists or not sig_exists:
        print('Either no Lennard-Jones epsilon or sigma found in'
              'in save filesys. Running OneDMin for params...')
        _run = True
    elif overwrite:
        print('User specified to overwrite parameters with new run...')
        _run = True
    else:
        _run = False

    if _run:

        # Obtain the geometry for the target and bath
        target_geo = get_geometry(
            tgt_info, thy_info,
            tgt_save_prefix,
            conf=conf)
        bath_geo = get_geometry(
            bath_info, thy_info,
            bath_save_prefix,
            conf=conf)

        # Run an instancw of 1DMin for each processor
        for idx in range(njobs):

            # Build run directory
            ljrunner.make_job_dirs(etrans_run_path, idx)

            # Write the input strings

            # Write the files

        # Submit the job
        print('\n\nRunning each OneDMin job...')
        ljrunner.submit_job(DRIVE_PATH, etrans_run_path, NJOBS)

    else:
        eps_path = etrans_save_fs[-1].file.epsilon.path(etrans_locs)
        sig_path = etrans_save_fs[-1].file.sigma.path(etrans_locs)
        print('- Lennard-Jones epsilon found at path {}'.format(eps_path))
        print('- Lennard-Jones sigma found at path {}'.format(sig_path))
        


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
