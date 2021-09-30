"""
Executes the automation part of 1DMin
"""

import autofile
import automol
import onedmin_io
import autorun
from mechanalyzer.inf import spc as sinfo
from mechanalyzer.inf import thy as tinfo
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter
from mechlib.amech_io._path import job_path


def run_onedmin(spc_name,
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
    _, tgt_cnf_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        spc_locs=tgt_info, thy_locs=tgt_mod_thy_info[1:])

    tgt_loc_info = filesys.mincnf.min_energy_conformer_locators(
        tgt_cnf_save_fs, tgt_mod_thy_info)
    _, tgt_cnf_save_path = tgt_loc_info

    # Build the bath conformer filesystem objects
    _, bath_cnf_save_fs = filesys.build_fs(
        run_prefix, save_prefix, 'CONFORMER',
        spc_locs=bath_info, thy_locs=bath_mod_thy_info[1:])

    # Build the target energy transfer filesystem objects
    etrans_save_fs = autofile.fs.energy_transfer(tgt_cnf_save_path)
    etrans_locs = bath_info + lj_mod_thy_info[1:4]

    # Calculate and save the Lennard-Jones parameters, if needed
    nsamp_needed = _nsamp_needed(
        etrans_save_fs, etrans_locs, etrans_keyword_dct)
    if nsamp_needed > 0:
        _runlj(nsamp_needed,
               lj_info, tgt_dct,
               lj_mod_thy_info, tgt_mod_thy_info, bath_mod_thy_info,
               tgt_cnf_save_fs, bath_cnf_save_fs,
               etrans_save_fs, etrans_locs,
               etrans_keyword_dct, run_prefix)
    else:
        epath = etrans_save_fs[-1].file.lennard_jones_epsilon.path(etrans_locs)
        spath = etrans_save_fs[-1].file.lennard_jones_sigma.path(etrans_locs)
        ioprinter.info_message(
            f'- Lennard-Jones epsilon found at path {epath}')
        ioprinter.info_message(
            f'- Lennard-Jones sigma found at path {spath}')


def _nsamp_needed(etrans_save_fs, etrans_locs, etrans_keyword_dct):
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
        nsamp_need = nsamp
    elif overwrite:
        ioprinter.info_message(
            'User specified to overwrite parameters with new run...')
        nsamp_need = nsamp
    else:
        inf_obj = etrans_save_fs[-1].file.info.read(etrans_locs)
        nsampd = inf_obj.nsamp
        if nsamp < nsampd:
            nsamp_need = nsampd - nsamp
        else:
            nsamp_need = 0

    return nsamp_need


def _runlj(nsamp_needed,
           lj_info, tgt_dct,
           lj_mod_thy_info, tgt_mod_thy_info, bath_mod_thy_info,
           tgt_cnf_save_fs, bath_cnf_save_fs,
           etrans_save_fs, etrans_locs,
           etrans_keyword_dct, run_prefix):
    """ Run the Lennard-Jones parameters
    """

    # Get property for target, bath, and combination
    tgt_geo = filesys.read.geometry(
        tgt_cnf_save_fs, tgt_mod_thy_info, conf=etrans_keyword_dct['conf'])
    bath_geo = filesys.read.geometry(
        bath_cnf_save_fs, bath_mod_thy_info, conf=etrans_keyword_dct['conf'])

    charge = sinfo.value(lj_info, 'charge')
    mult = sinfo.value(lj_info, 'mult')
    ich1, ich2 = sinfo.value(lj_info, 'inchi')
    fml_str = automol.formula.string(
        automol.formula.join(
            automol.inchi.formula(ich1),
            automol.inchi.formula(ich2)))

    # Set job parameters
    njobs = etrans_keyword_dct['njobs']
    nsamp_per_job = nsamp_needed // njobs + 1  # run slighty more than needed
    run_dir = job_path(run_prefix, 'ONEDMIN', 'LJ', fml_str)

    # Get the sp str from the submission scritpt (Remove bash?)
    sp_script_str = autorun.SCRIPT_DCT[tinfo.value(lj_mod_thy_info, 'program')]
    sp_script_str = ''.join(sp_script_str.splitlines()[1:])

    # Run OneDMin
    # prolly better to just get strings
    # need geoms geo_str = _output_str(jobdir, 'min_geoms.out')
    inp_strs, els_str, out_strs = autorun.onedmin.direct(
        sp_script_str, run_dir, nsamp_per_job, njobs,
        tgt_geo, bath_geo, lj_mod_thy_info, charge, mult,
        smin=tgt_dct['smin'], smax=tgt_dct['smax'], spin_method=1)

    # Parse out certain info
    epsilons, sigmas, geoms, ranseeds, version, input_str = _parse(
        inp_strs, out_strs)

    # Save information to filesys
    filesys.save.energy_transfer(
        etrans_save_fs, etrans_locs,
        epsilons, sigmas, geoms,
        ranseeds, version, input_str, els_str)


# Helpers
def _parse(input_strs_lst, output_strs_lst):
    """ read the program and version,
        maybe pull out the ranseeds for reproducibility?
    """

    # Get the main data from the output
    sigmas, epsilons, ranseeds, geos = (), (), (), ()
    for output_strs in output_strs_lst:
        sigmas, epsilons = onedmin_io.reader.lennard_jones(output_strs[1])
        ranseed = onedmin_io.reader.ranseed(output_strs[1])
        sigmas += sigmas
        epsilons += epsilons
        ranseeds += (ranseed,)

    version = onedmin_io.reader.program_version(output_strs_lst[1])

    # Get a random input str
    input_str = input_strs_lst[0]

    return epsilons, sigmas, geos, ranseeds, version, input_str
