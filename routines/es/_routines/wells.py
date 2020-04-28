"""
Find a TS from the grid as well as associated vdW wells
"""

import numpy
import automol
import elstruct
import autofile
from routines.es import _util as util
from routines.es import geom
from routines.es import runner as es_runner
from lib import filesys
from lib.phydat import phycon


def kick_from_saddle():
    """ Find the wells from kicking off the saddle point by
        changingthe  reaction coordinate some amount
        might have to be reaction class specific
    """


def find_with_irc():
    """ Try and use the wells by reading a computed irc, or computing an irc
        first and then reading it
    """


def find_vdw(ts_name, spc_dct, thy_info, ini_thy_info, vdw_params,
             nsamp_par, run_prefix, save_prefix,
             kickoff_size, kickoff_backward,
             overwrite):
    """ Find van der Waals structures for all the pairs of
        species in a reaction list.
        Fxn takes two species, performs a (random?) rotation,
        sticks them together and optimizes the combined geometry.
        Supposed to use the wells filesystem?
    """
    new_vdws = []
    _, opt_script_str, _, opt_kwargs = es_runner.par.run_qchem_par(
        *thy_info[:2])
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
        _, opt_script_str, _, opt_kwargs = es_runner.par.run_qchem_par(prog, method)

        geos = [(), ()]
        ntaudof = 0.
        for i, (nam, ich, chg, mul) in enumerate(zip(names, ichs, chgs, muls)):
            spc_info = [ich, chg, mul]
            orb_restr = filesys.inf.orbital_restriction(spc_info, ini_thy_info)
            ini_g = ini_thy_info[0:3]
            ini_g.append(orb_restr)
            orb_restr = filesys.inf.orbital_restriction(spc_info, thy_info)
            thy_info = thy_info[0:3]
            thy_info.append(orb_restr)
            spc_run_fs = autofile.fs.species(run_prefix)
            spc_run_fs[-1].create(spc_info)
            spc_run_path = spc_run_fs[-1].path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs[-1].create(spc_info)
            spc_save_path = spc_save_fs[-1].path(spc_info)

            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs[-1].create(thy_info[1:4])
            thy_run_path = thy_run_fs[-1].path(thy_info[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs[-1].create(thy_info[1:4])
            thy_save_path = thy_save_fs[-1].path(thy_info[1:4])
            run_fs = autofile.fs.run(thy_run_path)

            ini_thy_save_fs = autofile.fs.theory(spc_save_path)
            ini_thy_save_fs[-1].create(ini_thy_info[1:4])

            cnf_run_fs = autofile.fs.conformer(thy_run_path)
            cnf_save_fs = autofile.fs.conformer(thy_save_path)

            geo = geom.reference_geometry(
                spc_dct[nam], thy_info, ini_thy_info,
                thy_run_fs, thy_save_fs,
                ini_thy_save_fs,
                cnf_run_fs, cnf_save_fs,
                run_fs,
                kickoff_size=kickoff_size,
                kickoff_backward=kickoff_backward,
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
            spc_run_fs[-1].create(spc_info)
            spc_run_path = spc_run_fs[-1].path(spc_info)
            spc_save_fs = autofile.fs.species(save_prefix)
            spc_save_fs[-1].create(spc_info)
            spc_save_path = spc_save_fs[-1].path(spc_info)
            orb_restr = filesys.inf.orbital_restriction(spc_info, thy_info)
            thy_info = thy_info[0:3]
            thy_info.append(orb_restr)
            thy_run_fs = autofile.fs.theory(spc_run_path)
            thy_run_fs[-1].create(thy_info[1:4])
            thy_run_path = thy_run_fs[-1].path(thy_info[1:4])
            thy_save_fs = autofile.fs.theory(spc_save_path)
            thy_save_fs[-1].create(thy_info[1:4])
            thy_save_path = thy_save_fs[-1].path(thy_info[1:4])
            run_fs = autofile.fs.run(thy_run_path)

            # Generate reference geometry
            # Generate the z-matrix and sampling ranges
            es_runner.run_job(
                job=elstruct.Job.OPTIMIZATION,
                geom=geo,
                spc_info=spc_info,
                th_info=thy_info,
                run_fs=run_fs,
                script_str=opt_script_str,
                overwrite=overwrite,
                **opt_kwargs,
            )

            # Save info for the initial geometry (from ichi or fsave dir)
            ret = es_runner.read_job(
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
                thy_save_fs[-1].file.geometry.write(geo, thy_info[1:4])
                ene = elstruct.reader.energy(prog, method, out_str)
                if ene < min_ene:
                    min_ene = ene
                    print('ene test in vdw')
                    print(ene)
                    thy_save_fs[-1].file.energy.write(ene, thy_info[1:4])
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
                    cnf_save_fs[0].create()
                    cnf_run_fs[0].create()
                    tors_range_dct = {}
                    cinf_obj = autofile.system.info.conformer[0](
                        0, tors_range_dct)
                    cinf_obj.nsamp = 1
                    cnf_save_fs[0].file.info.write(cinf_obj)
                    locs_lst = cnf_save_fs[-1].existing()
                    if not locs_lst:
                        cid = autofile.system.generate_new_conformer_id()
                        locs = [cid]
                    else:
                        locs = locs_lst[0]
                    cnf_save_fs[-1].create(locs)
                    cnf_run_fs[-1].create(locs)
                    cnf_save_fs[-1].file.geometry_info.write(
                        inf_obj, locs)
                    cnf_save_fs[-1].file.geometry_input.write(
                        inp_str, locs)
                    cnf_save_fs[-1].file.energy.write(ene, locs)
                    cnf_save_fs[-1].file.geometry.write(geo, locs)
        if min_ene:
            new_vdws.append(vdw_name)

    return new_vdws


def fake_conf(thy_info, filesystem, inf=()):
    """ generate data to be used for a fake well I think?
    """
    cnf_save_fs = filesystem[5]
    cnf_run_fs = filesystem[4]
    thy_save_fs = filesystem[3]
    run_fs = filesystem[-1]
    thy_save_path = thy_save_fs[-1].path(thy_info[1:4])
    geo = thy_save_fs[-1].file.geometry.read(thy_info[1:4])
    if inf:
        inf_obj, ene = inf
    else:
        ene = thy_save_fs[-1].file.energy.read(thy_info[1:4])
        inf_obj = run_fs[0].file.info.read()
    tors_range_dct = {}
    cinf_obj = autofile.system.info.conformer_trunk(0, tors_range_dct)
    cinf_obj.nsamp = 1
    cnf_save_fs = autofile.fs.conformer(thy_save_path)
    cnf_save_fs[0].create()
    cnf_run_fs[0].create()
    cnf_save_fs[0].file.info.write(cinf_obj)
    cnf_run_fs[0].file.info.write(cinf_obj)
    locs_lst = cnf_save_fs[-1].existing()
    if not locs_lst:
        cid = autofile.system.generate_new_conformer_id()
        locs = [cid]
    else:
        locs = locs_lst[0]
    cnf_save_fs[-1].create(locs)
    cnf_run_fs[-1].create(locs)
    cnf_save_fs[-1].file.geometry_info.write(
        inf_obj, locs)
    cnf_run_fs[-1].file.geometry_info.write(
        inf_obj, locs)
    # method = inf_obj.method
    cnf_save_fs[-1].file.energy.write(ene, locs)
    cnf_run_fs[-1].file.energy.write(ene, locs)
    cnf_save_fs[-1].file.geometry.write(geo, locs)
    cnf_run_fs[-1].file.geometry.write(geo, locs)


def fake_geo_gen(tsk, thy_info, filesystem):
    """ generate data to be used for a fake well I think?
    """
    if 'conf' in tsk:
        fake_conf(thy_info, filesystem)
    if 'scan' in tsk:
        pass
    if 'tau' in tsk:
        pass
