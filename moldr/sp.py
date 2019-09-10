""" drivers for single point calculations
"""
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def run_single_point_energy(
        geo, spc_info, thy_level, sp_run_fs, sp_save_fs,
        script_str, overwrite, **kwargs):
    """ Find the energy for the minimum energy conformer
    """
    sp_run_fs.leaf.create(thy_level[1:4])
    sp_run_path = sp_run_fs.leaf.path(thy_level[1:4])

    sp_save_fs.leaf.create(thy_level[1:4])
    sp_save_path = sp_save_fs.leaf.path(thy_level[1:4])
    run_fs = autofile.fs.run(sp_run_path)

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


def run_minimum_energy_gradient(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, **kwargs):
    """ Find the gradient for the minimum energy conformer
    """

    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        run_fs = autofile.fs.run(min_cnf_run_path)
        print('Minimum energy conformer gradient')
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

            print(" - Reading gradient from output...")
            grad = elstruct.reader.gradient(inf_obj.prog, out_str)

            print(" - Saving gradient...")
            print(" - Save path: {}".format(min_cnf_save_path))
            cnf_save_fs.leaf.file.gradient_info.write(inf_obj, min_cnf_locs)
            cnf_save_fs.leaf.file.gradient_input.write(inp_str, min_cnf_locs)
            cnf_save_fs.leaf.file.gradient.write(grad, min_cnf_locs)


def run_minimum_energy_hessian(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, **kwargs):
    """ Find the hessian for the minimum energy conformer
    """

    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        run_fs = autofile.fs.run(min_cnf_run_path)
        moldr.driver.run_job(
            job='hessian',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs
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
                print('Freqs test')
                print(freqs)
                freqs = elstruct.util.harmonic_frequencies(
                    geo, hess, project=True)

                print(" - Saving hessian...")
                print(" - Save path: {}".format(min_cnf_save_path))
                cnf_save_fs.leaf.file.hessian_info.write(inf_obj, min_cnf_locs)
                cnf_save_fs.leaf.file.hessian_input.write(inp_str, min_cnf_locs)
                cnf_save_fs.leaf.file.hessian.write(hess, min_cnf_locs)
                cnf_save_fs.leaf.file.harmonic_frequencies.write(
                    freqs, min_cnf_locs)


def run_minimum_energy_vpt2(
        spc_info, thy_level, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, **kwargs):
    """  Run vpt2 for the minimum energy conformer
    """

    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
    if min_cnf_locs:
        min_cnf_run_path = cnf_run_fs.leaf.path(min_cnf_locs)
        min_cnf_save_path = cnf_save_fs.leaf.path(min_cnf_locs)
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        run_fs = autofile.fs.run(min_cnf_run_path)
        print('Minimum energy conformer vpt2')
        moldr.driver.run_job(
            job='vpt2',
            script_str=script_str,
            run_fs=run_fs,
            geom=geo,
            spc_info=spc_info,
            thy_level=thy_level,
            overwrite=overwrite,
            **kwargs
        )

        ret = moldr.driver.read_job(
            job='vpt2',
            run_fs=run_fs,
        )

        if ret is not None:
            inf_obj, inp_str, out_str = ret
            prog = inf_obj.prog

            vpt2_dict = elstruct.reader.vpt2(prog, out_str)
            anh_freqs = vpt2_dict['freqs']
            anh_zpe = vpt2_dict['zpe']
            xmat = vpt2_dict['xmat']
            vibrot_mat = vpt2_dict['vibrot']
            cent_dist = vpt2_dict['cent_dist']

            print(" - Reading anharmonic data from output...")
            # replace following with vpt2 information
#            print(" - Saving hessian...")
#            print(" - Save path: {}".format(min_cnf_save_path))
            cnf_save_fs.leaf.file.vpt2_info.write(inf_obj, min_cnf_locs)
            cnf_save_fs.leaf.file.vpt2_input.write(inp_str, min_cnf_locs)
            cnf_save_fs.leaf.file.anh_freqs.write(anh_freqs, min_cnf_locs)
#            cnf_save_fs.leaf.file.zpe.write(zpe, min_cnf_locs)
#            cnf_save_fs.leaf.file.harmonic_frequencies.write(freqs, min_cnf_locs)
            cnf_save_fs.leaf.file.anh_freqs.write(anh_freqs, min_cnf_locs)
            cnf_save_fs.leaf.file.xmat.write(xmat, min_cnf_locs)
#            cnf_save_fs.leaf.file.vibrot.write(vibrot, min_cnf_locs)
            cnf_save_fs.leaf.file.cent_dist.write(cent_dist, min_cnf_locs)


