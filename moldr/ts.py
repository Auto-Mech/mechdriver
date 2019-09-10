""" ts drivers
"""
from qcelemental import constants as qcc
import automol
import elstruct
import autofile
import moldr

WAVEN2KCAL = qcc.conversion_factor('wavenumber', 'kcal/mol')
EH2KCAL = qcc.conversion_factor('hartree', 'kcal/mol')

def ts_conformer_sampling(
        spc_info, geo, zma, tors_names, thy_level, thy_save_fs, cnf_run_fs, cnf_save_fs,
        script_str, overwrite, saddle=True, nsamp_par=(False, 3, 3, 1, 50, 50), 
        **kwargs):
    """ Find the minimum energy conformer by optimizing from nsamp random
    initial torsional states
    """
#    geo = ts_save_fs.trunk.file.geometry.read(thy_level[1:4])
#    zma = automol.geom.zmatrix(geo)
#    tors_names = automol.geom.zmatrix_torsion_coordinate_names(geo)
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(
        zma, tors_names)
    tors_range_dct = dict(zip(tors_names, tors_ranges))
    ich = spc_info[0]
    ntaudof = len(tors_names)
    nsamp = moldr.util.nsamp_init(nsamp_par, ntaudof)

    ts_save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
    )
#    print('zma test in conformer sampling')
#    print(automol.zmatrix.string(zma))

    moldr.conformer.run_conformers(
        zma=zma,
        spc_info=spc_info,
        thy_level=thy_level,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        **kwargs,
    )

    ts_save_conformers(
        cnf_run_fs=cnf_run_fs,
        cnf_save_fs=cnf_save_fs,
    )

    # save information about the minimum energy conformer in top directory
    min_cnf_locs = moldr.util.min_energy_conformer_locators(cnf_save_fs)
#    print('thy_save_fs:', thy_save_fs)
#    print('thy_level:', thy_level)
    if min_cnf_locs:
        geo = cnf_save_fs.leaf.file.geometry.read(min_cnf_locs)
        zma = cnf_save_fs.leaf.file.zmatrix.read(min_cnf_locs)

        thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
        thy_save_fs.leaf.file.zmatrix.write(zma, thy_level[1:4])


def ts_save_conformers(cnf_run_fs, cnf_save_fs):
    """ save the conformers that have been found so far
    """

    locs_lst = cnf_save_fs.leaf.existing()
    seen_geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                 for locs in locs_lst]
    seen_enes = [cnf_save_fs.leaf.file.energy.read(locs)
                 for locs in locs_lst]

    if not cnf_run_fs.trunk.exists():
        print("No conformers to save. Skipping...")
    else:
        for locs in cnf_run_fs.leaf.existing():
            run_path = cnf_run_fs.leaf.path(locs)
            run_fs = autofile.fs.run(run_path)
            print("Reading from conformer run at {}".format(run_path))

            ret = moldr.driver.read_job(job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)
                geo = elstruct.reader.opt_geometry(prog, out_str)
                zma = elstruct.reader.opt_zmatrix(prog, out_str)

#                zma = automol.geom.zmatrix(geo)
                unique = True

                for idx, geoi in enumerate(seen_geos):
                    enei = seen_enes[idx]
                    etol = 1.e-6
                    if automol.geom.almost_equal_coulomb_spectrum(
                            geo, geoi, rtol=1e-2):
                        if abs(ene-enei) < etol:
                            unique = False
                if not unique:
                    print(" - Geometry is not unique. Skipping...")
                else:
                    save_path = cnf_save_fs.leaf.path(locs)
                    print(" - Geometry is unique. Saving...")
                    print(" - Save path: {}".format(save_path))

                    cnf_save_fs.leaf.create(locs)
                    cnf_save_fs.leaf.file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs.leaf.file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs.leaf.file.energy.write(ene, locs)
                    cnf_save_fs.leaf.file.geometry.write(geo, locs)
                    cnf_save_fs.leaf.file.zmatrix.write(zma, locs)

                seen_geos.append(geo)
                seen_enes.append(ene)

        # update the conformer trajectory file
        locs_lst = cnf_save_fs.leaf.existing()
        if locs_lst:
            enes = [cnf_save_fs.leaf.file.energy.read(locs)
                    for locs in locs_lst]
            geos = [cnf_save_fs.leaf.file.geometry.read(locs)
                    for locs in locs_lst]

            traj = []
            for ene, geo in sorted(zip(enes, geos), key=lambda x: x[0]):
                comment = 'energy: {:>15.10f}'.format(ene)
                traj.append((comment, geo))

            traj_path = cnf_save_fs.trunk.file.trajectory.path()
            print("Updating conformer trajectory file at {}".format(traj_path))
            cnf_save_fs.trunk.file.trajectory.write(traj)


