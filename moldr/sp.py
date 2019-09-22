""" drivers for single point calculations
"""
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def run_energy(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Find the energy for the given structure
    """

    # prepare special file system since there may be many energies under same directory
    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    sp_run_fs = autofile.fs.single_point(geo_run_path)
    sp_save_fs = autofile.fs.single_point(geo_save_path)
    print('run_energy path test:', thy_level[1:4])
    sp_run_fs.leaf.create(thy_level[1:4])
    sp_run_path = sp_run_fs.leaf.path(thy_level[1:4])
    print('sp_run_path test:', sp_run_path)
    sp_save_fs.leaf.create(thy_level[1:4])
    sp_save_path = sp_save_fs.leaf.path(thy_level[1:4])
    print('sp_save_path test:', sp_run_path)
    run_fs = autofile.fs.run(sp_run_path)

    print('thy_level in run_energy:', thy_level) 
    moldr.driver.run_job(
        job='energy',
        script_str=script_str,
        run_fs=run_fs,
        geom=geo,
        spc_info=spc_info,
        thy_level=thy_level,
        overwrite=overwrite,
        **kwargs,
    )

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


def run_gradient(
        spc_info, thy_level, geo_run_fs, geo_save_fs, locs,
        script_str, overwrite, **kwargs):
    """ Determine the gradient for the geometry in the given location
    """

    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

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

    print('run hessian test:', locs)
    geo_run_path = geo_run_fs.leaf.path(locs)
    geo_save_path = geo_save_fs.leaf.path(locs)
    geo = geo_save_fs.leaf.file.geometry.read(locs)
    run_fs = autofile.fs.run(geo_run_path)

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
            anh = ()
        else:
            print(" - Reading hessian from output...")
            vpt2_dct = elstruct.reader.vpt2(inf_obj.prog, out_str)

            print(" - Saving hessian...")
            print(" - Save path: {}".format(geo_save_path))
            # geo_save_fs.leaf.file.vpt2_info.write(inf_obj, locs)
            geo_save_fs.leaf.file.vpt2_input.write(inp_str, locs)
            geo_save_fs.leaf.file.anharmonic_frequencies.write(vpt2_dct['freqs'], locs)
            geo_save_fs.leaf.file.anharmonic_zpve.write(vpt2_dct['zpve'], locs)
            geo_save_fs.leaf.file.anharmonicity_matrix.write(vpt2_dct['x_mat'], locs)
            geo_save_fs.leaf.file.vibro_rot_alpha_matrix.write(vpt2_dct['vibrot_mat'], locs)
            geo_save_fs.leaf.file.quartic_centrifugal_dist_consts.write(vpt2_dct['cent_dist_const'], locs)
