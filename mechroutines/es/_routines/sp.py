""" Routines for taking the geometry for a species or transition state
    conformer and calculating various molecular properties via some specified
    electronic structure method.A

    * Maybe put long description from run_energy here

    ** no reason to pass the geo and zma since the functions have the fs
       and locs and can therefore read them. the geo and zma are unneeded
"""

import automol
import elstruct
import autofile
import autorun
from phydat import phycon, symm
from mechlib.amech_io import printer as ioprinter
from mechroutines.es import runner as es_runner


# _JSON_SAVE = ['TAU']
_JSON_SAVE = []


def run_energy(zma, geo, spc_info, thy_info,
               geo_run_fs, geo_save_fs, locs,
               script_str, overwrite,
               retryfail=True, highspin=False, **kwargs):
    """ Assesses if an electronic energy exists in the CONFS/SP/THY layer
        of the save filesys for a species at the specified level of theory.
        If an energy does not exist, or if a user requests overwrite,
        the appropriate electronic struture calculation is set-up, launched,
        and then parsed within the run filesys, then that the energy and
        job input is written into the save filesys.

        :param zma: Z-Matrix for molecular structure to perform calculation
        :type zma: automol.zmat object
        :param geo: goemetry for molecular structure to perform calculation
        :type geo: automol.geom object
        :param spc_info:
        :type spc_info:
        :param thy_info:
        :type thy_info:
        :geo_run_fs: filesystem object for layer where structure is in RUN
        :geo_run_fs: autofile.fs object
        :geo_save_fs: filesystem object for layer where structure is in SAVE
        :geo_save_fs: autofile.fs object
        :param script_str: shell script for executing electronic structure job
        :type script_str: str
        :param overwrite:
    """

    geo_run_path = geo_run_fs[-1].path(locs)
    geo_save_path = geo_save_fs[-1].path(locs)

    # Prepare unique filesystem since many energies may be under same directory
    if not highspin:
        sp_run_fs = autofile.fs.single_point(geo_run_path)
        sp_save_fs = autofile.fs.single_point(geo_save_path)
    else:
        sp_run_fs = autofile.fs.high_spin(geo_run_path)
        sp_save_fs = autofile.fs.high_spin(geo_save_path)

    sp_run_path = sp_run_fs[-1].path(thy_info[1:4])
    sp_save_path = sp_save_fs[-1].path(thy_info[1:4])

    # Set input geom
    if geo is not None:
        job_geo = geo
    else:
        job_geo = zma

    exists = sp_save_fs[-1].file.energy.exists(thy_info[1:4])
    if not exists:
        ioprinter.info_message(
            'No energy found in save filesys. Running energy...')
        _run = True
    elif overwrite:
        ioprinter.info_message(
            'User specified to overwrite energy with new run...')
        _run = True
    else:
        _run = False

    if _run:

        # Add options matrix for energy runs for molpro
        if thy_info[0] == 'molpro2015':
            errs, optmat = es_runner.molpro_opts_mat(spc_info, geo)
        else:
            errs = ()
            optmat = ()

        sp_run_fs[-1].create(thy_info[1:4])
        run_fs = autofile.fs.run(sp_run_path)

        success, ret = es_runner.execute_job(
            job=elstruct.Job.ENERGY,
            script_str=script_str,
            run_fs=run_fs,
            geo=job_geo,
            spc_info=spc_info,
            thy_info=thy_info,
            errors=errs,
            options_mat=optmat,
            overwrite=overwrite,
            retryfail=retryfail,
            **kwargs,
        )

        if success:
            inf_obj, inp_str, out_str = ret

            ioprinter.info_message(" - Reading energy from output...")
            ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

            ioprinter.energy(ene)
            sp_save_fs[-1].create(thy_info[1:4])
            sp_save_fs[-1].file.input.write(inp_str, thy_info[1:4])
            sp_save_fs[-1].file.info.write(inf_obj, thy_info[1:4])
            sp_save_fs[-1].file.energy.write(ene, thy_info[1:4])
            ioprinter.save_energy(sp_save_path)

    else:
        ioprinter.existing_path('Energy', sp_save_path)
        ene = sp_save_fs[-1].file.energy.read(thy_info[1:4])
        ioprinter.energy(ene)


def run_gradient(zma, geo, spc_info, thy_info,
                 geo_run_fs, geo_save_fs, locs,
                 script_str, overwrite,
                 retryfail=True, **kwargs):
    """ Determine the gradient for the geometry in the given location
    """

    geo_run_path = geo_run_fs[-1].path(locs)
    geo_save_path = geo_save_fs[-1].path(locs)

    # Set input geom
    if geo is not None:
        job_geo = geo
    else:
        job_geo = automol.zmat.geometry(zma)
    is_atom = automol.geom.is_atom(job_geo)

    if not is_atom:

        if _json_database(geo_save_path):
            exists = geo_save_fs[-1].json.gradient.exists(locs)
        else:
            exists = geo_save_fs[-1].file.gradient.exists(locs)
        if not exists:
            ioprinter.info_message(
                'No gradient found in save filesys. Running gradient...')
            _run = True
        elif overwrite:
            ioprinter.info_message(
                'User specified to overwrite gradient with new run...')
            _run = True
        else:
            _run = False

        if _run:

            run_fs = autofile.fs.run(geo_run_path)

            success, ret = es_runner.execute_job(
                job=elstruct.Job.GRADIENT,
                script_str=script_str,
                run_fs=run_fs,
                geo=job_geo,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                retryfail=retryfail,
                **kwargs,
            )

            if success:
                inf_obj, inp_str, out_str = ret

                if is_atom:
                    grad = ()
                else:
                    ioprinter.info_message(
                        " - Reading gradient from output...")
                    grad = elstruct.reader.gradient(inf_obj.prog, out_str)

                    ioprinter.info_message(" - Saving gradient...")
                    if _json_database(geo_save_path):
                        geo_save_fs[-1].json.gradient_info.write(inf_obj, locs)
                        geo_save_fs[-1].json.gradient_input.write(
                            inp_str, locs)
                        geo_save_fs[-1].json.gradient.write(grad, locs)
                    else:
                        geo_save_fs[-1].file.gradient_info.write(inf_obj, locs)
                        geo_save_fs[-1].file.gradient_input.write(
                            inp_str, locs)
                        geo_save_fs[-1].file.gradient.write(grad, locs)
                    ioprinter.save_geo(geo_save_path)

        else:
            ioprinter.existing_path('Gradient', geo_save_path)

    else:
        ioprinter.info_message('Species is an atom. Skipping gradient task.')


def run_hessian(zma, geo, spc_info, thy_info,
                geo_run_fs, geo_save_fs, locs,
                script_str, overwrite,
                retryfail=True, **kwargs):
    """ Determine the hessian for the geometry in the given location
    """

    # Set the run filesystem information
    geo_run_path = geo_run_fs[-1].path(locs)
    geo_save_path = geo_save_fs[-1].path(locs)
    run_fs = autofile.fs.run(geo_run_path)

    # Set input geom
    # Warning using internal coordinates leads to inconsistencies with Gaussian
    # For this reason we only use Cartesians to generate the Hessian
    if geo is not None:
        job_geo = geo
    else:
        job_geo = automol.zmat.geometry(zma)
    is_atom = automol.geom.is_atom(job_geo)

    if not is_atom:

        if _json_database(geo_save_path):
            exists = geo_save_fs[-1].json.hessian.exists(locs)
        else:
            exists = geo_save_fs[-1].file.hessian.exists(locs)
        if not exists:
            ioprinter.info_message(
                'No Hessian found in save filesys. Running Hessian...')
            _run = True
        elif overwrite:
            ioprinter.info_message(
                'User specified to overwrite Hessian with new run...')
            _run = True
        else:
            _run = False

        if _run:

            run_fs = autofile.fs.run(geo_run_path)

            success, ret = es_runner.execute_job(
                job=elstruct.Job.HESSIAN,
                script_str=script_str,
                run_fs=run_fs,
                geo=job_geo,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                retryfail=retryfail,
                **kwargs,
            )

            if success:
                inf_obj, inp_str, out_str = ret

                ioprinter.info_message(" - Reading hessian from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)

                ioprinter.info_message(" - Saving Hessian...")
                if _json_database(geo_save_path):
                    geo_save_fs[-1].json.hessian_info.write(inf_obj, locs)
                    geo_save_fs[-1].json.hessian_input.write(inp_str, locs)
                    geo_save_fs[-1].json.hessian.write(hess, locs)
                else:
                    geo_save_fs[-1].file.hessian_info.write(inf_obj, locs)
                    geo_save_fs[-1].file.hessian_input.write(inp_str, locs)
                    geo_save_fs[-1].file.hessian.write(hess, locs)
                ioprinter.info_message(
                    f" - Save path: {geo_save_path}")

                if thy_info[0] == 'gaussian09':
                    _hess_grad(inf_obj.prog, out_str, geo_save_fs,
                               geo_save_path, locs, overwrite)
                _hess_freqs(geo, geo_save_fs,
                            geo_run_path, geo_save_path, locs, overwrite)

        else:
            ioprinter.existing_path('Hessian', geo_save_path)

            _hess_freqs(geo, geo_save_fs,
                        geo_run_path, geo_save_path, locs, overwrite)

    else:
        ioprinter.info_message('Species is an atom. Skipping Hessian task.')


def run_vpt2(zma, geo, spc_info, thy_info,
             geo_run_fs, geo_save_fs, locs,
             script_str, overwrite,
             retryfail=True, **kwargs):
    """ Perform vpt2 analysis for the geometry in the given location
    """

    # Set the run filesystem information
    geo_run_path = geo_run_fs[-1].path(locs)
    geo_save_path = geo_save_fs[-1].path(locs)
    run_fs = autofile.fs.run(geo_run_path)

    # Assess if symmetry needs to be broken for the calculation
    # Add porgram check because might only be issue for gaussian
    if spc_info[0] in symm.HIGH:
        if zma is not None:
            disp = symm.HIGH[spc_info[0]] * phycon.ANG2BOHR
            vals = automol.zmat.value_dictionary(zma)
            job_geo = automol.zmat.set_values_by_name(
                zma, {'R1': vals['R1'] + disp})
            is_atom = automol.geom.is_atom(
                automol.zmat.geometry(job_geo))
        else:
            ioprinter.warning_message(
                f'Need a zma for high-symmetry of {spc_info[0]}.',
                'Skipping task')
            job_geo = None
            is_atom = False
    else:
        if zma is not None:
            job_geo = zma
            is_atom = automol.geom.is_atom(automol.zmat.geometry(job_geo))
        else:
            job_geo = geo
            is_atom = automol.geom.is_atom(job_geo)

    if is_atom:
        ioprinter.info_message('Species is an atom, Skipping VPT2 task.')

    if job_geo is not None and not is_atom:

        exists = geo_save_fs[-1].file.anharmonicity_matrix.exists(locs)
        if not exists:
            ioprinter.info_message(
                'No X-Matrix found in save filesys. Running VPT2...')
            _run = True
        elif overwrite:
            ioprinter.info_message(
                'User specified to overwrite VPT2 info with new run...')
            _run = True
        else:
            _run = False

        if _run:

            success, ret = es_runner.execute_job(
                job=elstruct.Job.VPT2,
                script_str=script_str,
                run_fs=run_fs,
                geo=job_geo,
                spc_info=spc_info,
                thy_info=thy_info,
                overwrite=overwrite,
                retryfail=retryfail,
                **kwargs,
            )

            if success:
                inf_obj, inp_str, out_str = ret

                ioprinter.info_message(
                    " - Reading anharmonicities from output...")
                vpt2_dct = elstruct.reader.vpt2(inf_obj.prog, out_str)

                ioprinter.save_anharmonicity(geo_save_path)
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
                geo_save_fs[-1].file.cubic_force_constants.write(
                    vpt2_dct['cubic_fc'], locs)
                geo_save_fs[-1].file.quartic_force_constants.write(
                    vpt2_dct['quartic_fc'], locs)

        else:
            ioprinter.existing_path('VPT2 information', geo_save_path)


def run_prop(zma, geo, spc_info, thy_info,
             geo_run_fs, geo_save_fs, locs,
             script_str, overwrite,
             retryfail=True, **kwargs):
    """ Determine the properties in the given location
    """

    # Set input geom
    geo_run_path = geo_run_fs[-1].path(locs)
    geo_save_path = geo_save_fs[-1].path(locs)
    if geo is not None:
        job_geo = geo
    else:
        job_geo = automol.zmat.geometry(zma)

    if _json_database(geo_save_path):
        dmom_exists = geo_save_fs[-1].json.dipole_moment.exists(locs)
        polar_exists = geo_save_fs[-1].json.polarizability.exists(locs)
    else:
        dmom_exists = geo_save_fs[-1].file.dipole_moment.exists(locs)
        polar_exists = geo_save_fs[-1].file.polarizability.exists(locs)
    if not dmom_exists or not polar_exists:
        ioprinter.info_message(
            'Either no dipole moment or polarizability found in'
            'in save filesys. Running properties...')
        _run = True
    elif overwrite:
        ioprinter.info_message(
            'User specified to overwrite property with new run...')
        _run = True
    else:
        _run = False

    if _run:

        run_fs = autofile.fs.run(geo_run_path)

        success, ret = es_runner.execute_job(
            job=elstruct.Job.MOLPROP,
            script_str=script_str,
            run_fs=run_fs,
            geo=job_geo,
            spc_info=spc_info,
            thy_info=thy_info,
            overwrite=overwrite,
            retryfail=retryfail,
            **kwargs,
        )

        if success:
            inf_obj, _, out_str = ret

            ioprinter.info_message(" - Reading dipole moment from output...")
            dmom = elstruct.reader.dipole_moment(inf_obj.prog, out_str)
            ioprinter.info_message(" - Reading polarizability from output...")
            polar = elstruct.reader.polarizability(inf_obj.prog, out_str)

            ioprinter.debug_message('dip mom', dmom)
            ioprinter.debug_message('polar', polar)
            ioprinter.info_message(
                " - Saving dipole moment and polarizability...")
            if _json_database(geo_save_path):
                # geo_save_fs[-1].json.property_info.write(inf_obj, locs)
                # geo_save_fs[-1].json.property_input.write(
                #     inp_str, locs)
                geo_save_fs[-1].json.dipole_moment.write(dmom, locs)
                geo_save_fs[-1].json.polarizability.write(polar, locs)
            else:
                # geo_save_fs[-1].file.property_info.write(inf_obj, locs)
                # geo_save_fs[-1].file.property_input.write(
                #     inp_str, locs)
                geo_save_fs[-1].file.dipole_moment.write(dmom, locs)
                geo_save_fs[-1].file.polarizability.write(polar, locs)
            ioprinter.save_geo(geo_save_path)

    else:
        ioprinter.existing_path(
            'Dipole moment and polarizability',
            geo_save_path)


def _hess_freqs(geo, geo_save_fs, run_path, save_path, locs, overwrite):
    """ Calculate harmonic frequencies using Hessian
    """
    if _json_database(save_path):
        exists = geo_save_fs[-1].json.harmonic_frequencies.exists(locs)
    else:
        exists = geo_save_fs[-1].file.harmonic_frequencies.exists(locs)
    if not exists:
        ioprinter.info_message(
                'No harmonic frequencies found in save filesys...')
        _run = True
    elif overwrite:
        ioprinter.info_message(
                'User specified to overwrite frequencies with new run...')
        _run = True
    else:
        _run = False

    if _run:

        # Read the Hessian from the filesystem
        if _json_database(save_path):
            hess = geo_save_fs[-1].json.hessian.read(locs)
        else:
            hess = geo_save_fs[-1].file.hessian.read(locs)

        # Calculate and save the harmonic frequencies
        ioprinter.info_message(
                " - Calculating harmonic frequencies from Hessian...")
        script_str = autorun.SCRIPT_DCT['projrot']
        rt_freqs, _, rt_imags, _ = autorun.projrot.frequencies(
            script_str, run_path, [geo], [[]], [hess])
        rt_imags = tuple(-1 * imag_freq for imag_freq in rt_imags)
        freqs = sorted(rt_imags + rt_freqs)
        ioprinter.frequencies(freqs)
        ioprinter.geometry(geo)
        if _json_database(save_path):
            geo_save_fs[-1].json.harmonic_frequencies.write(freqs, locs)
        else:
            geo_save_fs[-1].file.harmonic_frequencies.write(freqs, locs)
        ioprinter.save_frequencies(save_path)

    else:
        ioprinter.existing_path('Harmonic frequencies', save_path)
        freqs = geo_save_fs[-1].file.harmonic_frequencies.read(locs)
        ioprinter.frequencies(freqs)


def _hess_grad(prog, out_str, geo_save_fs,
               save_path, locs, overwrite):
    """ Grab and save the gradient from a Hessian job if possible.
    """

    if _json_database(save_path):
        exists = geo_save_fs[-1].json.gradient.exists(locs)
    else:
        exists = geo_save_fs[-1].file.gradient.exists(locs)
    if not exists:
        ioprinter.info_message('No gradient found in save filesys...')
        _read = True
    elif overwrite:
        ioprinter.info_message(
            'User specified to overwrite gradient with new run...')
        _read = True
    else:
        _read = False

    if _read:

        # Read the Gradient from the electronic structure output
        ioprinter.info_message(
            " - Attempting to read gradient from Hessian from output...")
        grad = elstruct.reader.gradient(prog, out_str)

        if grad is not None:

            # Calculate and save the harmonic frequencies
            if _json_database(save_path):
                geo_save_fs[-1].json.gradient.write(grad, locs)
            else:
                geo_save_fs[-1].file.gradient.write(grad, locs)
            ioprinter.save_gradient(save_path)

        else:

            ioprinter.info_message(
                'No gradient found in Hessian job output to save.')

    else:
        ioprinter.existing_path('Gradient', save_path)


def _json_database(save_path):
    """Is this save path using a json style database (or directory style)
    """
    it_is = False
    for inst in _JSON_SAVE:
        if inst in save_path:
            it_is = True
    return it_is
