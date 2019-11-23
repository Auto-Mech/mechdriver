""" ts drivers
"""

import automol
import elstruct
import autofile
import moldr
import scripts.es


def reference_geometry(spcdct, thy_level, ini_thy_level,
                       geo_fs, ini_fs,
                       dist_info=(), overwrite=False):
    """ determine what to use as the reference geometry for all future runs
    If ini_thy_info refers to geometry dictionary then use that,
    otherwise values are from a hierarchy of:
    running level of theory, input level of theory, inchis.
    From the hierarchy an optimization is performed followed by a check for
    an imaginary frequency and then a conformer file system is set up.
    """

    thy_save_fs = geo_fs[3]
    ini_thy_save_fs = ini_fs[1]
    run_fs = geo_fs[-1]
    print('initializing geometry')
    geo = None
    geo_init = None
    # Check to see if geometry should be obtained from dictionary
    spc_info = [spcdct['ich'], spcdct['chg'], spcdct['mul']]
    if 'input_geom' in ini_thy_level:
        geom_obj = spcdct['geoobj']
        geo_init = geom_obj
        overwrite = True
        print('found initial geometry from geometry dictionary')
    else:
        # Check to see if geo already exists at running_theory
        if thy_save_fs.trunk.file.geometry.exists():
            thy_path = thy_save_fs.trunk.path()
            print(
                'getting reference geometry from {}'.format(thy_path))
            geo = thy_save_fs.trunk.file.geometry.read()
            zma = thy_save_fs.trunk.file.zmatrix.read()
            # print('geo:',automol.geom.string(geo))
        if not geo:
            if ini_thy_save_fs:
                if ini_thy_save_fs.trunk.file.geometry.exists():
                    thy_path = ini_thy_save_fs.trunk.path()
                    print(
                        'getting reference geometry from {}'.format(thy_path))
                    zma_init = ini_thy_save_fs.trunk.file.zmatrix.read()
                    geo_init = ini_thy_save_fs.trunk.file.geometry.read()
    if not geo and geo_init:
        _, opt_script_str, _, opt_kwargs = moldr.util.run_qchem_par(
            *thy_level[0:2])
        moldr.driver.run_job(
            job='optimization',
            script_str=opt_script_str,
            run_fs=run_fs,
            geom=zma_init,
            spc_info=spc_info,
            thy_level=thy_level,
            saddle=True,
            overwrite=overwrite,
            **opt_kwargs,
            )
        opt_ret = moldr.driver.read_job(
            job='optimization',
            run_fs=run_fs,
            )
        if opt_ret is not None:
            inf_obj, _, out_str = opt_ret
            prog = inf_obj.prog
            method = inf_obj.method
            ene = elstruct.reader.energy(prog, method, out_str)
            geo = elstruct.reader.opt_geometry(prog, out_str)
            zma = elstruct.reader.opt_zmatrix(prog, out_str)
        if geo:

            print(" - Saving...")
            print(" - Save path: {}".format(thy_save_fs.trunk.path()))

            thy_save_fs.trunk.file.energy.write(ene)
            thy_save_fs.trunk.file.geometry.write(geo)
            thy_save_fs.trunk.file.zmatrix.write(zma)

            scripts.es.run_single_conformer(
                spc_info, thy_level, geo_fs,
                overwrite, True, dist_info)
    return geo


def cas_options_1(spc_info, formula, num_act_elc, num_act_orb, high_mul):
    """ prepare cas options for multireference wavefunctions
    """

    elec_count = automol.formula.electron_count(formula)
    closed_orb = (elec_count-num_act_elc)//2
    occ_orb = closed_orb + num_act_orb
    closed_orb = (elec_count-num_act_elc)//2 - 2
    two_spin = spc_info[2]-1
    high_two_spin = high_mul - 1
    chg = spc_info[1]
    cas_opt = [
        elstruct.option.specify(
            elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(
            elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.WFN_, elec_count, 1, two_spin, chg)
        ]

    wfn_str = (
        "{{uhf,maxit=300;wf,{0},1,{1}}}\n"
        "{{multi,maxit=40;closed,{2};occ,{3};wf,{0},1,{4};orbprint,3}}"
    ).format(elec_count, high_two_spin, closed_orb, occ_orb, two_spin)

    return cas_opt, wfn_str


def cas_options_2(spc_info, formula, num_act_elc, num_act_orb, high_mul):
    """ prepare cas options for multireference wavefunctions
    """

    elec_count = automol.formula.electron_count(formula)
    closed_orb = (elec_count-num_act_elc)//2
    occ_orb = closed_orb + num_act_orb
    two_spin = spc_info[2]-1
    high_two_spin = high_mul - 1
    chg = spc_info[1]
    cas_opt = [
        elstruct.option.specify(
            elstruct.Option.Scf.MAXITER_, 40),
        elstruct.option.specify(
            elstruct.Option.Casscf.OCC_, occ_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.CLOSED_, closed_orb),
        elstruct.option.specify(
            elstruct.Option.Casscf.WFN_, elec_count, 1, two_spin, chg)
        ]

    wfn_str = (
        "{{uhf,maxit=300;wf,{0},1,{1}}}\n"
        "{{multi,maxit=40;closed,{2};occ,{3};wf,{0},1,{4};orbprint,3}}"
    ).format(elec_count, high_two_spin, closed_orb, occ_orb, two_spin)

    return cas_opt, wfn_str


def multiref_wavefunction_guess(high_mul, zma,
                                spc_info, thy_level,
                                casscf_options):
    """ prepare wavefunction template for multireference electronic structure calcs
    """

    charge = spc_info[1]
    mul = spc_info[2]
    basis = thy_level[2]
    prog = thy_level[0]
    prog = 'molpro2015'

    guess_str1 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=high_mul,
        method='hf',
        basis=basis,
        prog=prog,
        orb_restricted=False,
        mol_options=['nosym'],
        )
    guess_str1 += '\n\n'
    guess_str1 = '\n'.join(guess_str1.splitlines()[2:-6])

    guess_str2 = elstruct.writer.energy(
        geom=zma,
        charge=charge,
        mult=mul,
        method='casscf',
        basis=basis,
        prog=prog,
        orb_restricted=True,
        casscf_options=casscf_options[0],
        mol_options=['nosym'],
        )
    guess_str2 += '\n\n'
    guess_str2 = '\n'.join(guess_str2.splitlines()[2:-6])

    guess_str = guess_str1 + guess_str2 

    if len(casscf_options) > 1:

        guess_str3 = elstruct.writer.energy(
            geom=zma,
            charge=charge,
            mult=mul,
            method='casscf',
            basis=basis,
            prog=prog,
            orb_restricted=True,
            casscf_options=casscf_options[1],
            mol_options=['nosym'],
            )
        guess_str3 += '\n\n'
        guess_str3 = '\n'.join(guess_str3.splitlines()[2:])

        guess_str += guess_str3

    return guess_str


def _run_irc(
        ts_info, thy_level, ts_run_fs, ts_save_fs, locs,
        overwrite, new_grid=False, **opt_kwargs):
    """ Run the IRC
    """

    # Obtain saddle-point minimmum-energy conformer from filesystem
    ts_run_path = ts_run_fs.leaf.path(locs)
    # ts_save_path = ts_save_fs.leaf.path(locs)
    geo = ts_save_fs.leaf.file.geometry.read(locs)
    
    # Check if IRC run to desired specs
    # If Not run the IRC calculation
    for grid_idx, grid_val, run_prefix in zip(grid_idxs, grid_vals, run_prefixes):
        if not scn_save_fs.leaf.file.geometry.exists([['RX'], [grid_val]]) or overwrite:
            run_irc = True

    if run_irc:
        # Set up the IRC run filesystem
        run_fs = autofile.fs.run(ts_run_path)

        # Run the IRC in the forward and reverse direction
        for irc_direction in ('forward', 'reverse'):
            moldr.driver.run_job(
                job='irc',
                script_str=irc_script_str,
                run_fs=run_fs,
                geom=zma,
                spc_info=ts_info,
                thy_level=ref_level,
                overwrite=overwrite,
                irc_direction=irc_direction,
                **opt_kwargs,
                )

    # Read the IRC from the file system
    for irc_direction in ('forward', 'reverse'):
        opt_ret = moldr.driver.read_job(
            job=elstruct.Job.IRC,
            run_fs=run_fs,
        )

        if opt_ret is not None:
            inf_obj, _, out_str = opt_ret
            prog = inf_obj.prog
            geos, gras, hessians = elstruct.reader.irc_points(prog, out_str)
            enes = elstruct.reader.irc_energies(prog, out_str)
            coords = elstruct.reader.irc_coordinates(prog, out_str)

            print(" - Saving...")
            print(" - Save path: {}".format())

            dist_name = 'RX'
            for idx, coord in enumerate(coords):
                locs = [[dist_name], [coord]]
                # save_fs.leaf.file.energy.write(enes[idx], locs)
                # save_fs.leaf.file.geometry.write(geos[idx], locs)
                # save_fs.leaf.file.gradient.write(gras[idx], locs)
                # save_fs.leaf.file.hessian.write(hessians[idx], locs)
