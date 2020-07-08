"""
  Sampling routines
"""

import automol
import elstruct:


def torsional_sampling(zma, spc_info,
                       mod_thy_info, thy_save_fs,
                       cnf_run_fs, cnf_save_fs,
                       script_str, overwrite,
                       saddle=False, nsamp_par=(False, 3, 3, 1, 50, 50),
                       tors_names='',
                       two_stage=False, retryfail=True,
                       rxn_class='',
                       constraints=('tors', 'none'),
                       **kwargs):
    """ Run a torsional sampling procedure
    """

    tors_range_dct = _set_tors_range(tors_name_grps)
    nsamp = _set_nsamp(nsamp_par, ich, tors_range_dct, saddle)


def run_sampling(
        zma, spc_info, thy_info, nsamp, tors_range_dct,
        cnf_run_fs, cnf_save_fs, script_str, overwrite,
        saddle, two_stage, retryfail,
        **kwargs):
    """ run sampling algorithm to find conformers
    """

    cnf_save_fs[0].create()
    nsamp0 = nsamp
    inf_obj = autofile.schema.info_objects.conformer_trunk(0, tors_range_dct)
    nsampd = _find_nsampd(save_fs, run_fs)

    tot_samp = nsamp - nsampd
    print(' - Number of samples that have been currently run:', nsampd)
    print(' - Number of samples requested:', nsamp)

    if nsamp-nsampd > 0:
        print('\nRunning {} samples...'.format(nsamp-nsampd))
    while True:
        nsamp = nsamp0 - nsampd
        # Break the while loop if enough sampls completed
        if nsamp <= 0:
            print('Requested number of samples have been completed. '
                  'Conformer search complete.')
            break

        # Run the conformer sampling
        samp_zma = _generate_samp_zma(zma, tors_range_dct, nsampd)

        cid = autofile.schema.generate_new_conformer_id()
        locs = [cid]

        cnf_run_fs[-1].create(locs)
        cnf_run_path = cnf_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(cnf_run_path)

        print("Run {}/{}".format(nsamp, tot_samp))
        tors_names = list(tors_range_dct.keys())
        if two_stage and tors_names:
            print('Stage one beginning, holding the coordinates constant',
                  tors_names)
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                frozen_coordinates=[tors_names],
                saddle=saddle,
                retryfail=retryfail,
                **kwargs
            )
            print('Stage one success, reading for stage 2')
            ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                sinf_obj, _, out_str = ret
                prog = sinf_obj.prog
                samp_zma = elstruct.reader.opt_zmatrix(prog, out_str)
                print('Stage one success beginning stage two on', samp_zma)
                es_runner.run_job(
                    job=elstruct.Job.OPTIMIZATION,
                    script_str=script_str,
                    run_fs=run_fs,
                    geom=samp_zma,
                    spc_info=spc_info,
                    thy_info=thy_info,
                    overwrite=overwrite,
                    saddle=saddle,
                    retryfail=False,
                    **kwargs
                )
        else:
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                script_str=script_str,
                run_fs=run_fs,
                geom=samp_zma,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                saddle=saddle,
                retryfail=retryfail,
                **kwargs
            )

        # Determine nsampd
        nsampd = _find_nsampd(save_fs, run_fs)
        nsampd += 1

        # Update nsampd values in filesys
        inf_obj.nsamp = nsampd
        save_fs[0].file.info.write(inf_obj)
        run_fs[0].file.info.write(inf_obj)


def _generate_samp_zma(zma, tors_range_dct, nsampd):
    """ Generate a sample ZMA
    """
    if nsampd > 0:
        samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
    else:
        samp_zma = zma

    return samp_zma


def _set_tors_range(tors_name_grps):
    """ Set the range for all of the torsions
    """
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(tors_name_grps)
    tors_range_dct = dict(zip(
        tuple(grp[0] for grp in tors_name_grps), tors_ranges))

    return tors_range_dct


# Determine variables for the number of samples
def _set_nsamp(nsamp_par, ich, tors_range_dct, saddle):
    """ Determine the number of samples to run
    """

    if tors_range_dct:

        # Generate nsamp using num_tors and parameters
        if not saddle:
            gra = automol.inchi.graph(ich)
            ntaudof = len(
                automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        else:
            ntaudof = len((tors_range_dct.keys()))

        # Set nsamp
        if nsamp_par[0]:
            nsamp = min(nsamp_par[1] + nsamp_par[2] * nsamp_par[3]**ntaudof,
                        nsamp_par[4])
        else:
            nsamp = nsamp_par[5]

    else:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    return nsamp


def _find_nsampd(save_fs, run_fs):
    """ Determine the number of coordinates that have been sampled
        using info files from the filesystem
    """

    if save_fs[0].file.info.exists():
        inf_obj_s = save_fs[0].file.info.read()
        nsampd = inf_obj_s.nsamp
    elif run_fs[0].file.info.exists():
        inf_obj_r = run_fs[0].file.info.read()
        nsampd = inf_obj_r.nsamp
    else:
        nsampd = 0

    return nsampd
