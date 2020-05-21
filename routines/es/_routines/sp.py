""" es_runners for single point calculations
"""

import automol
import elstruct
import autofile
from routines.es import runner as es_runner
from lib.phydat import phycon
from lib.phydat import symm


def run_energy(zma, geo, spc_info, thy_info,
               geo_save_fs, geo_run_path, geo_save_path, locs,
               script_str, overwrite,
               retryfail=True, **kwargs):
    """ Find the energy for the given structure
    """

    # geo_save_fs and locs unneeded for this
    _, _ = geo_save_fs, locs

    # Prepare unique filesystem since many energies may be under same directory
    sp_run_fs = autofile.fs.single_point(geo_run_path)
    sp_save_fs = autofile.fs.single_point(geo_save_path)
    sp_run_fs[-1].create(thy_info[1:4])
    sp_run_path = sp_run_fs[-1].path(thy_info[1:4])
    sp_save_fs[-1].create(thy_info[1:4])
    run_fs = autofile.fs.run(sp_run_path)

    # Set input geom
    if geo is not None:
        job_geo = geo
    else:
        job_geo = zma

    if not sp_save_fs[-1].file.energy.exists(thy_info[1:4]) or overwrite:

        print('No energy found in save filesys. Running energy...')
        # Add options matrix for energy runs for molpro
        if thy_info[0] == 'molpro2015':
            errs, optmat = es_runner.par.set_molpro_options_mat(spc_info, geo)
        else:
            errs = ()
            optmat = ()

        es_runner.run_job(
            job='energy',
            script_str=script_str,
            run_fs=run_fs,
            geom=job_geo,
            spc_info=spc_info,
            thy_info=thy_info,
            errors=errs,
            options_mat=optmat,
            overwrite=overwrite,
            retryfail=retryfail,
            **kwargs,
        )

        ret = es_runner.read_job(
            job='energy',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading energy from output...")
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving energy...")
            sp_save_fs[-1].file.input.write(inp_str, thy_info[1:4])
            sp_save_fs[-1].file.info.write(inf_obj, thy_info[1:4])
            sp_save_fs[-1].file.energy.write(ene, thy_info[1:4])

    else:
        print('Energy found and saved previously at {}'.format(
            sp_save_fs[-1].file.energy.path(thy_info[1:4])))


def run_gradient(zma, geo, spc_info, thy_info,
                 geo_save_fs, geo_run_path, geo_save_path, locs,
                 script_str, overwrite,
                 retryfail=True, **kwargs):
    """ Determine the gradient for the geometry in the given location
    """

    # Set input geom
    if geo is not None:
        job_geo = geo
    else:
        job_geo = automol.zmatrix.geometry(zma)
    is_atom = automol.geom.is_atom(job_geo)

    if not is_atom:

        if not geo_save_fs[-1].file.gradient.exists(locs) or overwrite:

            run_fs = autofile.fs.run(geo_run_path)

            print('No gradient found in save filesys. Running gradient...')
            es_runner.run_job(
                job='gradient',
                script_str=script_str,
                run_fs=run_fs,
                geom=job_geo,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                retryfail=retryfail,
                **kwargs,
            )

            ret = es_runner.read_job(
                job='gradient',
                run_fs=run_fs,
            )

            if ret is not None:
                inf_obj, inp_str, out_str = ret

                if is_atom:
                    grad = ()
                else:
                    print(" - Reading gradient from output...")
                    grad = elstruct.reader.gradient(inf_obj.prog, out_str)

                    print(" - Saving gradient...")
                    print(" - Save path: {}".format(geo_save_path))
                    geo_save_fs[-1].file.gradient_info.write(inf_obj, locs)
                    geo_save_fs[-1].file.gradient_input.write(inp_str, locs)
                    geo_save_fs[-1].file.gradient.write(grad, locs)

        else:
            print('Gradient found and saved previously at {}'.format(
                geo_save_fs[-1].file.gradient.path(locs)))

    else:
        print('Species is an atom, Skipping Gradient task.')


def run_hessian(zma, geo, spc_info, thy_info,
                geo_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite,
                retryfail=True, **kwargs):
    """ Determine the hessian for the geometry in the given location
    """

    # Set the run filesystem information
    run_fs = autofile.fs.run(geo_run_path)

    # if prog == 'molpro2015':
    #     geo = hess_geometry(out_str)
    #     scn_save_fs[-1].file.geometry.write(geo, locs)

    # Set input geom
    # Warning using internal coordinates leads to inconsistencies with Gaussian
    # For this reason we only use Cartesians to generate the Hessian
    if geo is not None:
        job_geo = geo
    else:
        job_geo = automol.zmatrix.geometry(zma)
    is_atom = automol.geom.is_atom(job_geo)

    if not is_atom:

        if not geo_save_fs[-1].file.hessian.exists(locs) or overwrite:

            run_fs = autofile.fs.run(geo_run_path)

            print('No Hessian found in save filesys. Running Hessian...')
            es_runner.run_job(
                job='hessian',
                script_str=script_str,
                run_fs=run_fs,
                geom=job_geo,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                retryfail=retryfail,
                **kwargs,
            )

            ret = es_runner.read_job(
                job='hessian',
                run_fs=run_fs,
            )

            if ret is not None:
                inf_obj, inp_str, out_str = ret

                print(" - Reading hessian from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=False)

                print(" - Saving Hessian...")
                print(" - Save path: {}".format(geo_save_path))
                geo_save_fs[-1].file.hessian_info.write(inf_obj, locs)
                geo_save_fs[-1].file.hessian_input.write(inp_str, locs)
                geo_save_fs[-1].file.hessian.write(hess, locs)
                geo_save_fs[-1].file.harmonic_frequencies.write(freqs, locs)

        else:
            print('Hessian found and saved previously at {}'.format(
                geo_save_fs[-1].file.hessian.path(locs)))

    else:
        print('Species is an atom, Skipping Hessian task.')


def run_vpt2(zma, geo, spc_info, thy_info,
             geo_save_fs, geo_run_path, geo_save_path, locs,
             script_str, overwrite,
             retryfail=True, **kwargs):
    
    """ Perform vpt2 analysis for the geometry in the given location
    """

    # Set the run filesystem information
    run_fs = autofile.fs.run(geo_run_path)

    # Assess if symmetry needs to be broken for the calculation
    # Add porgram check because might only be issue for gaussian
    if spc_info[0] in symm.HIGH:
        if zma is not None:
            disp = symm.HIGH[spc_info[0]] * phycon.ANG2BOHR
            vals = automol.zmatrix.values(zma)
            job_geo = automol.zmatrix.set_values(
                zma, {'R1': vals['R1'] + disp})
            is_atom = automol.geom.is_atom(
                automol.zmatrix.geometry(job_geo))
        else:
            print('Need a zma for high-symmetry of {}.'.format(spc_info[0]),
                  'Skipping task')
            job_geo = None
            is_atom = False
    else:
        if zma is not None:
            job_geo = zma
            is_atom = automol.geom.is_atom(automol.zmatrix.geometry(job_geo))
        else:
            job_geo = geo
            is_atom = automol.geom.is_atom(job_geo)

    if is_atom:
        print('Species is an atom, Skipping VPT2 task.')

    if job_geo is not None and not is_atom:

        run_vpt2_job = bool(
            not geo_save_fs[-1].file.anharmonicity_matrix.exists(locs) or
            overwrite)
        if run_vpt2_job:
            print('Running vpt2')
            es_runner.run_job(
                job='vpt2',
                script_str=script_str,
                run_fs=run_fs,
                geom=job_geo,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                retryfail=retryfail,
                **kwargs,
            )

            ret = es_runner.read_job(
                job='vpt2',
                run_fs=run_fs,
            )

            if ret is not None:
                inf_obj, inp_str, out_str = ret

                # if not geo_save_fs[-1].file.hessian.exists(locs):
                #     print(" - No Hessian in filesys.",
                #           "Reading it from output...")
                #     hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                #     print(" - Saving Hessian...")
                #     print(" - Save path: {}".format(geo_save_path))
                #     geo_save_fs[-1].file.hessian_info.write(inf_obj, locs)
                #     geo_save_fs[-1].file.hessian_input.write(inp_str, locs)
                #     geo_save_fs[-1].file.hessian.write(hess, locs)

                print(" - Reading anharmonicities from output...")
                vpt2_dct = elstruct.reader.vpt2(inf_obj.prog, out_str)
                # hess = elstruct.reader.hessian(inf_obj.prog, out_str)

                # Write the VPT2 file specifying the Fermi Treatments
                # fermi_treatment = '{} Defaults'.format(inf_obj.prog)
                # vpt2_inf_obj = autofile.schema.info_objects.vpt2(
                #     fermi_treatment=fermi_treatment)

                print(" - Saving anharmonicities...")
                print(" - Save path: {}".format(geo_save_path))
                # geo_save_fs[-1].file.vpt2_info.write(inf_obj, locs)
                geo_save_fs[-1].file.vpt2_input.write(inp_str, locs)
                geo_save_fs[-1].file.anharmonic_frequencies.write(
                    vpt2_dct['freqs'], locs)
                geo_save_fs[-1].file.anharmonic_zpve.write(
                    vpt2_dct['zpve'], locs)
                geo_save_fs[-1].file.vibro_rot_alpha_matrix.write(
                    vpt2_dct['vibrot_mat'], locs)
                geo_save_fs[-1].file.quartic_centrifugal_dist_consts.write(
                    vpt2_dct['cent_dist_const'], locs)
                geo_save_fs[-1].file.anharmonicity_matrix.write(
                    vpt2_dct['x_mat'], locs)

        else:
            print('VPT2 information found and saved previously at {}'.format(
                geo_save_fs[-1].file.anharmonicity_matrix.exists(locs)))
