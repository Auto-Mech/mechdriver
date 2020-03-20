""" drivers for single point calculations
"""

import automol
import elstruct
import autofile

# New libs
from lib.runner import par as runpar
from lib.runner import driver
from lib.phydat import symm, phycon

def run_energy(geo, spc_info, thy_level,
               geo_save_fs, geo_run_path, geo_save_path, locs,
               script_str, overwrite, **kwargs):
    """ Find the energy for the given structure
    """

    # Prepare unique filesystem since many energies may be under same directory
    sp_run_fs = autofile.fs.single_point(geo_run_path)
    sp_save_fs = autofile.fs.single_point(geo_save_path)
    sp_run_fs[-1].create(thy_level[1:4])
    sp_run_path = sp_run_fs[-1].path(thy_level[1:4])
    sp_save_fs[-1].create(thy_level[1:4])
    run_fs = autofile.fs.run(sp_run_path)

    if not sp_save_fs[-1].file.energy.exists(thy_level[1:4]) or overwrite:

        # Add options matrix for energy runs for molpro
        if thy_level[0] == 'molpro2015':
            errors, options_mat = runpar.set_molpro_options_mat(spc_info, geo)
        else:
            errors = ()
            options_mat = ()

        driver.run_job(
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

    ret = driver.read_job(
        job='energy',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        print(" - Reading energy from output...")
        ene = elstruct.reader.energy(inf_obj.prog, inf_obj.method, out_str)

        print(" - Saving energy...")
        sp_save_fs[-1].file.input.write(inp_str, thy_level[1:4])
        sp_save_fs[-1].file.info.write(inf_obj, thy_level[1:4])
        sp_save_fs[-1].file.energy.write(ene, thy_level[1:4])


def run_gradient(geo, spc_info, thy_level,
                 geo_save_fs, geo_run_path, geo_save_path, locs,
                 script_str, overwrite, **kwargs):
    """ Determine the gradient for the geometry in the given location
    """

    run_fs = autofile.fs.run(geo_run_path)
    if not geo_save_fs[-1].file.gradient.exists(locs) or overwrite:
        print('Running gradient')
        driver.run_job(
            job='gradient',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

    ret = driver.read_job(
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
            geo_save_fs[-1].file.gradient_info.write(inf_obj, locs)
            geo_save_fs[-1].file.gradient_input.write(inp_str, locs)
            geo_save_fs[-1].file.gradient.write(grad, locs)


def run_hessian(geo, spc_info, thy_level,
                geo_save_fs, geo_run_path, geo_save_path, locs,
                script_str, overwrite, **kwargs):
    """ Determine the hessian for the geometry in the given location
    """
    # if prog == 'molpro2015':
    #     geo = hess_geometry(out_str)
    #     scn_save_fs[-1].file.geometry.write(geo, locs)

    run_fs = autofile.fs.run(geo_run_path)
    if not geo_save_fs[-1].file.hessian.exists(locs) or overwrite:
        print('Running hessian')
        driver.run_job(
            job='hessian',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

    ret = driver.read_job(
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

            print(" - Saving Hessian...")
            print(" - Save path: {}".format(geo_save_path))
            geo_save_fs[-1].file.hessian_info.write(inf_obj, locs)
            geo_save_fs[-1].file.hessian_input.write(inp_str, locs)
            geo_save_fs[-1].file.hessian.write(hess, locs)
            geo_save_fs[-1].file.harmonic_frequencies.write(freqs, locs)


def run_vpt2(zma, spc_info, thy_level,
             geo_save_fs, geo_run_path, geo_save_path, locs,
             script_str, overwrite, **kwargs):
    """ Perform vpt2 analysis for the geometry in the given location
    """

    # Set the filesystem information
    run_fs = autofile.fs.run(geo_run_path)

    # Assess if symmetry needs to be broken for the calculation
    # Add porgram check because might only be issue for gaussian
    if spc_info[0] in symm.HIGH:
        disp = symm.HIGH[spc_info[0]] * phycon.ANG2BOHR
        vals = automol.zmatrix.values(zma)
        zma = automol.zmatrix.set_values(zma, {'R1': vals['R1'] + disp})

    run_vpt2 = bool(
        not geo_save_fs[-1].file.anharmonicity_matrix.exists(locs) or 
        not automol.geom.is_atom(geo) or
        overwrite)
    if run_vpt2:
        print('Running vpt2')
        driver.run_job(
            job='vpt2',
            script_str=script_str,
            run_fs=run_fs,
            geom=zma,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs,
        )

    ret = driver.read_job(
        job='vpt2',
        run_fs=run_fs,
    )

    if ret is not None:
        inf_obj, inp_str, out_str = ret

        if automol.geom.is_atom(geo):
            pass
        else:

            if not geo_save_fs[-1].file.hessian.exists(locs):
                print(" - No Hessian in filesys. Reading it from output...")
                hess = elstruct.reader.hessian(inf_obj.prog, out_str)
                print(" - Saving Hessian...")
                print(" - Save path: {}".format(geo_save_path))
                geo_save_fs[-1].file.hessian_info.write(inf_obj, locs)
                geo_save_fs[-1].file.hessian_input.write(inp_str, locs)
                geo_save_fs[-1].file.hessian.write(hess, locs)

            print(" - Reading anharmonicities from output...")
            vpt2_dct = elstruct.reader.vpt2(inf_obj.prog, out_str)
            hess = elstruct.reader.hessian(inf_obj.prog, out_str)

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
