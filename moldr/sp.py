""" drivers for single point calculations
"""
import automol
import elstruct
import autofile
import moldr


def run_energy(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Find the energy for the given structure
    """

    # Prepare unique filesystem since many energies may be under same directory
    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    sp_run_fs = autofile.fs.single_point(geo_run_path)
    sp_save_fs = autofile.fs.single_point(geo_save_path)
    sp_run_fs.leaf.create(thy_level[1:4])
    sp_run_path = sp_run_fs.leaf.path(thy_level[1:4])
    sp_save_fs.leaf.create(thy_level[1:4])
    sp_save_path = sp_save_fs.leaf.path(thy_level[1:4])
    run_fs = autofile.fs.run(sp_run_path)

    ret = moldr.driver.read_job(
        job='energy',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        print(" - Reading energy from output...")
        ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

        print(" - Saving energy...")
        print(" - Save path: {}".format(sp_save_path))
        sp_save_fs.leaf.file.energy.write(ene, thy_level[1:4])
        sp_save_fs.leaf.file.input.write(inp_str, thy_level[1:4])
        sp_save_fs.leaf.file.info.write(inf_obj, thy_level[1:4])

    if not sp_save_fs.leaf.file.energy.exists(thy_level[1:4]) or overwrite:

        # Add options matrix for energy runs for molpro
        if thy_level[0] == 'molpro2015':
            errors, options_mat = _set_molpro_options_mat(spc_info, geo)
        else:
            errors = ()
            options_mat = ()

        moldr.driver.run_job(
            job='energy',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            errors=errors,
            options_mat=options_mat,
            overwrite=overwrite,
            **kwargs,
        )

        ret = moldr.driver.read_job(
            job='energy',
            run_fs=run_fs,
        )
        print('ret in molpro options_mat:', ret)

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            print(" - Reading energy from output...")
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            print(" - Saving energy...")
            sp_save_fs.leaf.file.input.write(inp_str, thy_level[1:4])
            sp_save_fs.leaf.file.info.write(inf_obj, thy_level[1:4])
            sp_save_fs.leaf.file.energy.write(ene, thy_level[1:4])


def _set_molpro_options_mat(spc_info, geo):
    """ prepare the errors and options mat to perform successive
        single-point energy calculations in Molpro when the RHF fails to
        converge. This currently only works for doublets.
    """

    # Get the nelectrons, spins, and orbitals for the wf card
    formula = automol.geom.formula(geo)
    elec_count = automol.formula.electron_count(formula)
    num_act_elc = 1
    num_act_orb = 1
    closed_orb = (elec_count - num_act_elc) // 2
    occ_orb = closed_orb + num_act_orb
    two_spin = spc_info[2] - 1

    # Build the strings UHF and CASSCF wf card and set the errors and options
    uhf_str = (
        "{{uhf,maxit=300;wf,{0},1,{1};orbprint,3}}"
    ).format(elec_count, two_spin)
    cas_str = (
        "{{casscf,maxit=40;"
        "closed,{0};occ,{1};wf,{2},1,{3};canonical;orbprint,3}}"
    ).format(closed_orb, occ_orb, elec_count, two_spin)

    errors = [elstruct.Error.SCF_NOCONV]
    options_mat = [
        [{'gen_lines': {2: [uhf_str]}},
         {'gen_lines': {2: [cas_str]}},
         {'gen_lines': {2: [cas_str]}}
         ]
    ]

    return errors, options_mat


def run_gradient(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Determine the gradient for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

    ret = moldr.driver.read_job(
        job='gradient',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        if automol.geom.is_atom(geo):
            grad = ()
        else:
            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(geo_save_path))
            geo_save_fs.leaf.file.gradient_info.write(inf_obj, locs)
            geo_save_fs.leaf.file.gradient_input.write(inp_str, locs)
            geo_save_fs.leaf.file.gradient.write(grad, locs)

    if not geo_save_fs.leaf.file.gradient.exists(locs) or overwrite:
        print('Running gradient')
        moldr.driver.run_job(
            job='gradient',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = moldr.driver.read_job(
            job='gradient',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            if automol.geom.is_atom(geo):
                grad = ()
            else:
                print(" - Reading gradient from output...")
                grad = elstruct.reader.gradient(inf_obj.prog, out_str)

                print(" - Saving gradient...")
                print(" - Save path: {}".format(geo_save_path))
                geo_save_fs.leaf.file.gradient_info.write(inf_obj, locs)
                geo_save_fs.leaf.file.gradient_input.write(inp_str, locs)
                geo_save_fs.leaf.file.gradient.write(grad, locs)


def run_hessian(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Determine the hessian for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

    ret = moldr.driver.read_job(
        job='hessian',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        if automol.geom.is_atom(geo):
            freqs = ()
        else:
            print(" - Reading hessian from output...")
            hess = elstruct.reader.hessian(inf_obj.prog, out_str)
            freqs = elstruct.util.harmonic_frequencies(
                geo, hess, project=False)

            print(" - Saving hessian...")
            print(" - Save path: {}".format(geo_save_path))
            geo_save_fs.leaf.file.hessian_info.write(inf_obj, locs)
            geo_save_fs.leaf.file.hessian_input.write(inp_str, locs)
            geo_save_fs.leaf.file.hessian.write(hess, locs)
            geo_save_fs.leaf.file.harmonic_frequencies.write(freqs, locs)

    if not geo_save_fs.leaf.file.hessian.exists(locs) or overwrite:
        print('Running hessian')
        moldr.driver.run_job(
            job='hessian',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

        ret = moldr.driver.read_job(
            job='hessian',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            if automol.geom.is_atom(geo):
                freqs = ()
            else:
                print(" - Reading hessian from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=False)

                print(" - Saving hessian...")
                print(" - Save path: {}".format(geo_save_path))
                geo_save_fs.leaf.file.hessian_info.write(inf_obj, locs)
                geo_save_fs.leaf.file.hessian_input.write(inp_str, locs)
                geo_save_fs.leaf.file.hessian.write(hess, locs)
                geo_save_fs.leaf.file.harmonic_frequencies.write(freqs, locs)


def run_vpt2(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Perform vpt2 analysis for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

    ret = moldr.driver.read_job(
        job='vpt2',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        if automol.geom.is_atom(geo):
            pass
        else:
            print(" - Reading anharmonicities from output...")
            vpt2_dct = elstruct.reader.vpt2(inf_obj.prog, out_str)

            print(" - Saving anharmonicities...")
            print(" - Save path: {}".format(geo_save_path))
            # geo_save_fs.leaf.file.vpt2_info.write(inf_obj, locs)
            geo_save_fs.leaf.file.vpt2_input.write(inp_str, locs)
            geo_save_fs.leaf.file.anharmonic_frequencies.write(
                vpt2_dct['freqs'], locs)
            geo_save_fs.leaf.file.anharmonic_zpve.write(
                vpt2_dct['zpve'], locs)
            geo_save_fs.leaf.file.anharmonicity_matrix.write(
                vpt2_dct['x_mat'], locs)
            geo_save_fs.leaf.file.vibro_rot_alpha_matrix.write(
                vpt2_dct['vibrot_mat'], locs)
            geo_save_fs.leaf.file.quartic_centrifugal_dist_consts.write(
                vpt2_dct['cent_dist_const'], locs)

    if not geo_save_fs.leaf.file.anharmnicity_matrix.exists(locs) or overwrite:

        print('Running vpt2')
        moldr.driver.run_job(
            job='vpt2',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )
        ret = moldr.driver.read_job(
            job='vpt2',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret

            if automol.geom.is_atom(geo):
                pass
            else:
                print(" - Reading anharmonicities from output...")
                vpt2_dct = elstruct.reader.vpt2(inf_obj.prog, out_str)

                print(" - Saving anharmonicities...")
                print(" - Save path: {}".format(geo_save_path))
                # geo_save_fs.leaf.file.vpt2_info.write(inf_obj, locs)
                geo_save_fs.leaf.file.vpt2_input.write(inp_str, locs)
                geo_save_fs.leaf.file.anharmonic_frequencies.write(
                    vpt2_dct['freqs'], locs)
                geo_save_fs.leaf.file.anharmonic_zpve.write(
                    vpt2_dct['zpve'], locs)
                geo_save_fs.leaf.file.anharmonicity_matrix.write(
                    vpt2_dct['x_mat'], locs)
                geo_save_fs.leaf.file.vibro_rot_alpha_matrix.write(
                    vpt2_dct['vibrot_mat'], locs)
                geo_save_fs.leaf.file.quartic_centrifugal_dist_consts.write(
                    vpt2_dct['cent_dist_const'], locs)
