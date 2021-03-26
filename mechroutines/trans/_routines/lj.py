"""
Executes the automation part of 1DMin
"""

import statistics
import autofile
from autorun import run_script
from mechroutines.trans._routines import _geom as geom
from mechroutines.trans._routines import _gather as gather
from mechroutines.trans.runner import lj as lj_runner
from mechlib import filesys
from mechlib import amech_io
from mechlib.amech_io import printer as ioprinter


def onedmin(spc_name,
            spc_dct, thy_dct, etrans_keyword_dct,
            run_prefix, save_prefix):
    """ Run the task
    """

    bath_name = etrans_keyword_dct['bath']

    tgt_dct, bath_dct = spc_dct[spc_name], spc_dct[bath_name]
    tgt_info = filesys.inf.get_spc_info(tgt_dct)
    bath_info = filesys.inf.get_spc_info(bath_dct)
    lj_info = filesys.inf.combine_spc_info(tgt_info, bath_info)

    # Build the modified thy objs
    inp_thy_info = filesys.inf.get_es_info(
        etrans_keyword_dct['inplvl'], thy_dct)
    run_thy_info = filesys.inf.get_es_info(
        etrans_keyword_dct['runlvl'], thy_dct)
    tgt_mod_thy_info = filesys.inf.modify_orb_restrict(
        tgt_info, inp_thy_info)
    bath_mod_thy_info = filesys.inf.modify_orb_restrict(
        bath_info, inp_thy_info)
    lj_mod_thy_info = filesys.inf.modify_orb_restrict(
        lj_info, run_thy_info)

    # Build the target conformer filesystem objects
    tgt_cnf_run_fs, tgt_cnf_save_fs = build_fs(
        run_prefix, save_prefix, 'RING_CONFORMER',
        spc_locs=tgt_info, thy_locs=tgt_mod_thy_info[1:])

    tgt_loc_info = filesys.mincnf.min_energy_ring_conformer_locators(
        tgt_cnf_save_fs, tgt_mod_thy_info)
    locs, paths = tgt_loc_info
    tgt_min_rng_locs, tgt_min_cnf_locs = locs
    tgt_min_rng_path, tgt_min_cnf_path = paths

    # Create run fs if that directory has been deleted to run the jobs
    tgt_rng_run_fs[-1].create(tgt_min_rng_locs)
    tgt_rng_run_path = tgt_rng_run_fs[-1].path(tgt_min_rng_locs)
    tgt_cnf_run_fs, _ = build_fs(
       tgti_rng_run_path, tgt_rng_run_path, 'CONFORMER')
    tgt_cnf_run_fs[-1].create(tgt_min_cnf_locs)
    tgt_cnf_run_path = tgt_cnf_run_fs[-1].path(tgt_min_cnf_locs)
    tgt_cnf_save_path = tgt_min_cnf_path

    # Get options from the dct or es options lst
    # tgt_cnf_run_fs[-1].create(tgt_min_cnf_locs)
    # tgt_cnf_run_path = filesys.build.cnf_paths_from_locs(
    #     tgt_cnf_run_fs, [tgt_min_cnf_locs])[0]





    # Build the target energy transfer filesystem objects
    etrans_run_fs = autofile.fs.energy_transfer(tgt_cnf_run_path)
    etrans_save_fs = autofile.fs.energy_transfer(tgt_cnf_save_path)
    etrans_locs = bath_info + lj_mod_thy_info[1:4]

    # Build the bath conformer filesystem objects
    # _, bath_thy_save_path = filesys.build.spc_thy_fs_from_root(
    #     save_prefix, bath_info, bath_mod_thy_info)
    # ioprinter.debug_message('bath path', bath_thy_save_path)
    # bath_cnf_save_fs = autofile.fs.conformer(bath_thy_save_path)

    # Calculate and save the Lennard-Jones parameters, if needed
    run_needed, nsamp_needed = _need_run(
        etrans_save_fs, etrans_locs, etrans_keyword_dct)
    if run_needed:
        _runlj(nsamp_needed,
               lj_info, lj_mod_thy_info,
               tgt_mod_thy_info, bath_mod_thy_info,
               tgt_cnf_save_fs, bath_cnf_save_fs,
               etrans_run_fs, etrans_locs,
               etrans_keyword_dct)
        _savelj(etrans_run_fs, etrans_save_fs, etrans_locs,
                etrans_keyword_dct)
    else:
        epath = etrans_save_fs[-1].file.lennard_jones_epsilon.path(etrans_locs)
        spath = etrans_save_fs[-1].file.lennard_jones_sigma.path(etrans_locs)
        ioprinter.info_message(
            '- Lennard-Jones epsilon found at path {}'.format(epath))
        ioprinter.info_message(
            '- Lennard-Jones sigma found at path {}'.format(spath))


def _need_run(etrans_save_fs, etrans_locs, etrans_keyword_dct):
    """ Check if job needs to run
    """

    nsamp = etrans_keyword_dct['nsamp']
    overwrite = etrans_keyword_dct['overwrite']

    ex1 = etrans_save_fs[-1].file.lennard_jones_epsilon.exists(etrans_locs)
    ex2 = etrans_save_fs[-1].file.lennard_jones_sigma.exists(etrans_locs)
    if not ex1 or not ex2:
        ioprinter.info_message(
            'Either no Lennard-Jones epsilon or sigma found in',
            'save filesys. Running OneDMin for params...')
        run = True
        nsamp_need = nsamp
    elif overwrite:
        ioprinter.info_message(
            'User specified to overwrite parameters with new run...')
        run = True
        nsamp_need = nsamp
    else:
        inf_obj = etrans_save_fs[-1].file.info.read(etrans_locs)
        nsampd = inf_obj.nsamp
        if nsamp < nsampd:
            run = True
            nsamp_need = nsampd - nsamp
        else:
            run = False
            nsamp_need = 0

    return run, nsamp_need


def _runlj(nsamp_needed,
           lj_info, lj_mod_thy_info,
           tgt_mod_thy_info, bath_mod_thy_info,
           tgt_cnf_save_fs, bath_cnf_save_fs,
           etrans_run_fs, etrans_locs,
           etrans_keyword_dct):
    """ Run the Lennard-Jones parameters
    """

    # Pull stuff from dct
    njobs = etrans_keyword_dct['njobs']
    smin = etrans_keyword_dct['smin']
    smax = etrans_keyword_dct['smax']
    conf = etrans_keyword_dct['conf']

    # Determine the number of samples per job
    nsamp_per_job = nsamp_needed // njobs

    # Set the path to the executable
    onedmin_exe_path = '/lcrc/project/CMRP/amech/OneDMin/build'

    # Obtain the geometry for the target and bath
    tgt_geo = geom.get_geometry(
        tgt_cnf_save_fs, tgt_mod_thy_info, conf=conf)
    bath_geo = geom.get_geometry(
        bath_cnf_save_fs, bath_mod_thy_info, conf=conf)

    # Set the path to the etrans lead fs
    etrans_run_path = etrans_run_fs[-1].path(etrans_locs)

    # Build the run directory
    onedmin_run_path = lj_runner.build_rundir(etrans_run_path)

    # Run an instancw of 1DMin for each processor
    for idx in range(njobs):

        # Build run directory
        onedmin_job_path = lj_runner.make_jobdir(onedmin_run_path, idx)

        # Write the input files
        xyz1_str, xyz2_str = lj_runner.write_xyz(tgt_geo, bath_geo)

        elstruct_inp_str, elstruct_sub_str = lj_runner.write_elstruct_inp(
            lj_info, lj_mod_thy_info)

        onedmin_str = lj_runner.write_input(
            nsamp_per_job, smin=smin, smax=smax,
            target_name='target.xyz', bath_name='bath.xyz')

        input_strs = (
            xyz1_str, xyz2_str,
            elstruct_inp_str, elstruct_sub_str,
            onedmin_str)
        input_names = (
            'target.xyz', 'bath.xyz',
            'qc.mol', 'ene.x',
            'input.dat')
        inp = tuple(zip(input_strs, input_names))
        amech_io.writer.write_files(
            inp, onedmin_job_path, exe_names=('ene.x'))

    # Write the batch submission script for each instance
    onedmin_sub_str = lj_runner.write_onedmin_sub(
        njobs, onedmin_run_path, onedmin_exe_path,
        exe_name='onedmin-dd-molpro.x')
    sub_inp = ((onedmin_sub_str, 'build.sh'),)
    amech_io.writer.write_files(
        sub_inp, onedmin_run_path, exe_names=('build.sh'))

    # Submit the all of the OneDMin jobs
    ioprinter.info_message('Running each OneDMin job...', newline=2)
    run_script(onedmin_sub_str, onedmin_run_path)


def _savelj(etrans_run_fs, etrans_save_fs, etrans_locs,
            etrans_keyword_dct):
    """ Save the Lennard-Jones parameters
    """

    # Read the dictionary
    ljpotential = etrans_keyword_dct['pot']

    # Set the run path to read the files
    etrans_run_path = etrans_run_fs[-1].path(etrans_locs)

    # Read any epsilons and sigma currently in the filesystem
    ioprinter.info_message(
        'Reading Lennard-Jones parameters and Geoms from filesystem...',
        newline=1)
    fs_geoms, fs_epsilons, fs_sigmas = gather.read_filesys(
        etrans_save_fs, etrans_locs)
    gather.print_lj_parms(fs_sigmas, fs_epsilons)

    # Read the lj from all the output files
    ioprinter.info_message(
        'Reading Lennard-Jones parameters and Geoms from output...',
        newline=1)
    run_geoms, run_epsilons, run_sigmas = gather.read_output(etrans_run_path)
    gather.print_lj_parms(run_sigmas, run_epsilons)

    # Read the program and version for onedmin
    prog_version = gather.prog_version(etrans_run_path)

    # Add the lists from the two together
    geoms = fs_geoms + run_geoms
    sigmas = fs_sigmas + run_sigmas
    epsilons = fs_epsilons + run_epsilons

    # Average the sigma and epsilon values
    if geoms and sigmas and epsilons:

        assert len(geoms) == len(sigmas) == len(epsilons), (
            'Number of geoms, sigmas, and epsilons not the same'
        )

        avg_sigma = statistics.mean(sigmas)
        avg_epsilon = statistics.mean(epsilons)
        nsampd = len(sigmas)
        ioprinter.info_message(
            'Average Sigma to save [unit]:', avg_sigma, newline=1)
        ioprinter.info_message('Average Epsilont to save [unit]:', avg_epsilon)
        ioprinter.info_message('Number of values = ', nsampd)

        # Update the trajectory file
        traj = []
        for geo, eps, sig in zip(geoms, epsilons, sigmas):
            comment = 'Epsilon: {}   Sigma: {}'.format(eps, sig)
            traj.append((comment, geo))

        # Write the info obj
        inf_obj = autofile.schema.info_objects.lennard_jones(
            nsampd, potential=ljpotential,
            program='OneDMin', version=prog_version)

        # Set up the electronic structure input file
        onedmin_inp_str = '<ONEDMIN INP>'
        els_inp_str = '<ELSTRUCT INP>'

        # Write the params to the save file system
        etrans_save_fs[-1].file.lj_input.write(onedmin_inp_str, etrans_locs)
        etrans_save_fs[-1].file.info.write(inf_obj, etrans_locs)
        etrans_save_fs[-1].file.molpro_inp_file.write(els_inp_str, etrans_locs)
        etrans_save_fs[-1].file.epsilon.write(avg_epsilon, etrans_locs)
        etrans_save_fs[-1].file.sigma.write(avg_sigma, etrans_locs)
        etrans_save_fs[1].file.trajectory.write(traj, etrans_locs)
