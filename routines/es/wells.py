"""
Find a TS from the grid as well as associated vdW wells
"""

import numpy
import automol
import elstruct
import autofile
from routines.es import util
from routines.es import geom
from lib.runner import runpar
from lib.runner import driver
from lib.runner import script
from lib.filesystem import orb as fsorb
from lib.phydat import phycon


def find_vdw(ts_name, spc_dct, thy_info, ini_thy_info, vdw_params,
             nsamp_par, run_prefix, save_prefix,
             kickoff_size, kickoff_backward,
             overwrite):
    """ Find van der Waals structures for all the pairs of
        species in a reaction list
    """
    projrot_script_str = script.PROJROT
    new_vdws = []
    _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(*thy_info[:2])
    mul = spc_dct[ts_name]['low_mul']
    vdw_names_lst = []
    if vdw_params[0]:
        vdw_names_lst.append([sorted(spc_dct[ts_name]['reacs']), mul, 'r'])
    if vdw_params[1]:
        vdw_names_lst.append([sorted(spc_dct[ts_name]['prods']), mul, 'p'])

    for names, ts_mul, label in vdw_names_lst:
        if len(names) < 2:
            print('Cannot find van der Waals well for unimolecular',
                  'reactant or product')
        ichs = list(map(lambda name: spc_dct[name]['ich'], names))
        chgs = list(map(lambda name: spc_dct[name]['chg'], names))
        muls = list(map(lambda name: spc_dct[name]['mul'], names))

        # theory
        prog = thy_info[0]
        method = thy_info[1]
        _, opt_script_str, _, opt_kwargs = runpar.run_qchem_par(prog, method)

        geos = [(), ()]
        ntaudof = 0.
        for i, (nam, ich, chg, mul) in enumerate(zip(names, ichs, chgs, muls)):
            spc_info = [ich, chg, mul]
            orb_restr = fsorb.orbital_restriction(spc_info, ini_thy_info)
            ini_thy_level = ini_thy_info[0:3]
            ini_thy_level.append(orb_restr)
            orb_restr = fsorb.orbital_restriction(spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info)
            spc_run_path = spc_run_fs.leaf.path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
            run_fs = autofile.fs.run(thy_run_path)

            ini_thy_save_fs = autofile.fs.theory(spc_save_path)
            ini_thy_save_fs.leaf.create(ini_thy_level[1:4])

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)

            ini_filesys = [None, ini_thy_save_fs]
            filesys = [spc_run_fs, spc_save_fs, thy_run_fs, thy_save_fs,
                       cnf_run_fs, cnf_save_fs, None, None,
                       None, None, run_fs]
            geo = geom.reference_geometry(
                spc_dct[nam], thy_level, ini_thy_level,
                filesys, ini_filesys,
                kickoff_size=kickoff_size,
                kickoff_backward=kickoff_backward,
                projrot_script_str=projrot_script_str,
                overwrite=overwrite)
            geos[i] = geo
            gra = automol.geom.graph(geo)
            ntaudof += len(
                automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
        nsamp = util.nsamp_init(nsamp_par, ntaudof)
        geo1, geo2 = geos
        geo1 = automol.geom.mass_centered(geo1)
        geo2 = automol.geom.mass_centered(geo2)
        min_ene = 0.
        for idx in range(int(nsamp)):
            print('Optimizing vdw geometry {}/{}'.format(idx+1, nsamp))
            angs1 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs2 = numpy.multiply(
                numpy.random.rand(3), [1*numpy.pi, 2*numpy.pi, 2*numpy.pi])
            angs12 = numpy.multiply(
                numpy.random.rand(2), [1*numpy.pi, 2*numpy.pi])
            geo1 = automol.geom.euler_rotated(geo1, *angs1)
            geo2 = automol.geom.euler_rotated(geo2, *angs2)
            dist_cutoff = 3.0 * phycon.ANG2BOHR

            geo = automol.geom.join(geo1, geo2, dist_cutoff, *angs12)
            print("Species: {}".format('+'.join(names)))
            print('vdw starting geometry')
            print(automol.geom.xyz_string(geo))

            # Set up the filesystem
            ich = automol.inchi.recalculate(automol.inchi.join(ichs))
            chg = sum(chgs)
            mul = ts_mul
            spc_info = (ich, chg, mul)
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs.leaf.create(spc_info)
            spc_run_path = spc_run_fs.leaf.path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs.leaf.create(spc_info)
            spc_save_path = spc_save_fs.leaf.path(spc_info)
            orb_restr = fsorb.orbital_restriction(spc_info, thy_info)
            thy_level = thy_info[0:3]
            thy_level.append(orb_restr)
            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs.leaf.create(thy_level[1:4])
            thy_run_path = thy_run_fs.leaf.path(thy_level[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs.leaf.create(thy_level[1:4])
            thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
            run_fs = autofile.fs.run(thy_run_path)

            # Generate reference geometry
            # Generate the z-matrix and sampling ranges
            driver.run_job(
                job=elstruct.Job.OPTIMIZATION,
                geom=geo,
                spc_info=spc_info,
                thy_level=thy_level,
                run_fs=run_fs,
                script_str=opt_script_str,
                overwrite=overwrite,
                **opt_kwargs,
            )

            # Save info for the initial geometry (from ichi or fsave dir)
            ret = driver.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                print('Saving reference geometry')
                print(" - Save path: {}".format(thy_save_path))

                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                geo = elstruct.reader.opt_geometry(prog, out_str)
                print('vdw ending geometry')
                print(automol.geom.xyz_string(geo))
                thy_save_fs.leaf.file.geometry.write(geo, thy_level[1:4])
                ene = elstruct.reader.energy(prog, method, out_str)
                if ene < min_ene:
                    min_ene = ene
                    print('ene test in vdw')
                    print(ene)
                    thy_save_fs.leaf.file.energy.write(ene, thy_level[1:4])
                    print('Saving reference geometry')
                    print(" - Save path: {}".format(thy_save_path))
                    vdw_name = label + ts_name.replace('ts', 'vdw')
                    spc_dct[vdw_name] = spc_dct[ts_name].copy()
                    spc_dct[vdw_name]['ich'] = ich
                    spc_dct[vdw_name]['mul'] = mul
                    spc_dct[vdw_name]['chg'] = chg
                    spc_dct[vdw_name]['dist_info'][1] = dist_cutoff

                    # Make a fake conformer
                    cnf_save_fs = autofile.fs.conformer(thy_save_path)
                    cnf_run_fs = autofile.fs.conformer(thy_run_path)
                    cnf_save_fs.trunk.create()
                    cnf_run_fs.trunk.create()
                    tors_range_dct = {}
                    cinf_obj = autofile.system.info.conformer_trunk(
                        0, tors_range_dct)
                    cinf_obj.nsamp = 1
                    cnf_save_fs.trunk.file.info.write(cinf_obj)
                    locs_lst = cnf_save_fs.leaf.existing()
                    if not locs_lst:
                        cid = autofile.system.generate_new_conformer_id()
                        locs = [cid]
                    else:
                        locs = locs_lst[0]
                    cnf_save_fs.leaf.create(locs)
                    cnf_run_fs.leaf.create(locs)
                    cnf_save_fs.leaf.file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs.leaf.file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs.leaf.file.energy.write(ene, locs)
                    cnf_save_fs.leaf.file.geometry.write(geo, locs)
        if min_ene:
            new_vdws.append(vdw_name)

    return new_vdws


def fake_conf(thy_level, filesys, inf=()):
    """ generate data to be used for a fake well I think?
    """
    cnf_save_fs = filesys[5]
    cnf_run_fs = filesys[4]
    thy_save_fs = filesys[3]
    run_fs = filesys[-1]
    thy_save_path = thy_save_fs.leaf.path(thy_level[1:4])
    geo = thy_save_fs.leaf.file.geometry.read(thy_level[1:4])
    if inf:
        inf_obj, ene = inf
    else:
        ene = thy_save_fs.leaf.file.energy.read(thy_level[1:4])
        inf_obj = run_fs.trunk.file.info.read()
    tors_range_dct = {}
    cinf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
    cinf_obj.nsamp = 1
    cnf_save_fs = autofile.fs.conformer(thy_save_path)
    cnf_save_fs.trunk.create()
    cnf_run_fs.trunk.create()
    cnf_save_fs.trunk.file.info.write(cinf_obj)
    cnf_run_fs.trunk.file.info.write(cinf_obj)
    locs_lst = cnf_save_fs.leaf.existing()
    if not locs_lst:
        cid = autofile.system.generate_new_conformer_id()
        locs = [cid]
    else:
        locs = locs_lst[0]
    cnf_save_fs.leaf.create(locs)
    cnf_run_fs.leaf.create(locs)
    cnf_save_fs.leaf.file.geometry_info.write(
        inf_obj, locs)
    cnf_run_fs.leaf.file.geometry_info.write(
        inf_obj, locs)
    # method = inf_obj.method
    cnf_save_fs.leaf.file.energy.write(ene, locs)
    cnf_run_fs.leaf.file.energy.write(ene, locs)
    cnf_save_fs.leaf.file.geometry.write(geo, locs)
    cnf_run_fs.leaf.file.geometry.write(geo, locs)


def fake_geo_gen(tsk, thy_level, filesys):
    """ generate data to be used for a fake well I think?
    """
    if 'conf' in tsk:
        fake_conf(thy_level, filesys)
    if 'scan' in tsk:
        pass
    if 'tau' in tsk:
        pass
