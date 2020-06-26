""" es_runners
"""

# import itertools
import numpy
import automol
import elstruct
import autofile
from routines.es._routines import _util as util
from routines.es import runner as es_runner
from lib import filesys
from lib.phydat import phycon


def tau_sampling(zma, ref_ene, spc_info, tors_name_grps, nsamp_par,
                 mod_thy_info,
                 tau_run_fs, tau_save_fs,
                 script_str, overwrite,
                 saddle=False, **opt_kwargs):
    """ Sample over torsions optimizing all other coordinates
    """

    # Read the geometry from the initial filesystem and set sampling
    tors_ranges = automol.zmatrix.torsional_sampling_ranges(tors_name_grps)
    tors_range_dct = dict(zip(
        tuple(grp[0] for grp in tors_name_grps), tors_ranges))
    gra = automol.zmatrix.graph(zma)
    ntaudof = len(automol.graph.rotational_bond_keys(gra, with_h_rotors=False))
    nsamp = util.nsamp_init(nsamp_par, ntaudof)

    # Run through tau sampling process
    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        mod_thy_info=mod_thy_info
    )

    run_tau(
        zma=zma,
        spc_info=spc_info,
        thy_info=mod_thy_info,
        nsamp=nsamp,
        tors_range_dct=tors_range_dct,
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        script_str=script_str,
        overwrite=overwrite,
        saddle=saddle,
        **opt_kwargs,
    )

    save_tau(
        tau_run_fs=tau_run_fs,
        tau_save_fs=tau_save_fs,
        mod_thy_info=mod_thy_info
    )

    print('Assessing the convergence of the Monte Carlo Partition Function...')
    assess_pf_convergence(tau_save_fs, ref_ene)


def run_tau(zma, spc_info, thy_info, nsamp, tors_range_dct,
            tau_run_fs, tau_save_fs, script_str, overwrite,
            saddle, **kwargs):
    """ run sampling algorithm to find tau dependent geometries
    """
    if not tors_range_dct:
        print("No torsional coordinates. Setting nsamp to 1.")
        nsamp = 1

    tau_save_fs[0].create()

    vma = automol.zmatrix.var_(zma)
    if tau_save_fs[0].file.vmatrix.exists():
        existing_vma = tau_save_fs[0].file.vmatrix.read()
        assert vma == existing_vma
    tau_save_fs[0].file.vmatrix.write(vma)
    idx = 0
    nsamp0 = nsamp
    inf_obj = autofile.schema.info_objects.tau_trunk(0, tors_range_dct)
    while True:
        if tau_save_fs[0].file.info.exists():
            inf_obj_s = tau_save_fs[0].file.info.read()
            nsampd = inf_obj_s.nsamp
        elif tau_run_fs[0].file.info.exists():
            inf_obj_r = tau_run_fs[0].file.info.read()
            nsampd = inf_obj_r.nsamp
        else:
            nsampd = 0

        nsamp = nsamp0 - nsampd
        if nsamp <= 0:
            print('Reached requested number of samples. '
                  'Tau sampling complete.')
            break

        print("    New nsamp is {:d}.".format(nsamp))

        samp_zma, = automol.zmatrix.samples(zma, 1, tors_range_dct)
        tid = autofile.schema.generate_new_tau_id()
        locs = [tid]

        tau_run_fs[-1].create(locs)
        tau_run_prefix = tau_run_fs[-1].path(locs)
        run_fs = autofile.fs.run(tau_run_prefix)

        idx += 1
        print("Run {}/{}".format(idx, nsamp0))

        print('\nChecking if ZMA has high repulsion...')
        if _low_repulsion_struct(zma, samp_zma):
            print('ZMA fine.')
            # es_runner.run_job(
                # job=elstruct.Job.OPTIMIZATION,
                # script_str=script_str,
                # run_fs=run_fs,
                # geom=samp_zma,
                # spc_info=spc_info,
                # thy_info=thy_info,
                # saddle=saddle,
                # overwrite=overwrite,
                # frozen_coordinates=tors_range_dct.keys(),
                # **kwargs
            # )
        else:
            print('repulsive ZMA:')
            inp_str = elstruct.writer.optimization(
                geom=samp_zma,
                charge=spc_info[1],
                mult=spc_info[2],
                method=thy_info[1],
                basis=thy_info[2],
                prog=thy_info[0],
                orb_type=thy_info[3],
                mol_options=['nosym'],
                frozen_coordinates=tors_range_dct.keys(),
            )
            tau_run_fs[-1].file.geometry_input.write(inp_str, locs)
            print('geometry for bad ZMA at',tau_run_fs[-1].path(locs))

        nsampd += 1
        inf_obj.nsamp = nsampd
        tau_save_fs[0].file.info.write(inf_obj)
        tau_run_fs[0].file.info.write(inf_obj)


def save_tau(tau_run_fs, tau_save_fs, mod_thy_info):
    """ save the tau dependent geometries that have been found so far
    """

    saved_geos = [tau_save_fs[-1].file.geometry.read(locs)
                  for locs in tau_save_fs[-1].existing()]

    if not tau_run_fs[0].exists():
        print("No tau geometries to save. Skipping...")
    else:
        for locs in tau_run_fs[-1].existing():
            run_path = tau_run_fs[-1].path(locs)
            run_fs = autofile.fs.run(run_path)

            print("Reading from tau run at {}".format(run_path))

            ret = es_runner.read_job(
                job=elstruct.Job.OPTIMIZATION, run_fs=run_fs)
            if ret:
                inf_obj, inp_str, out_str = ret
                prog = inf_obj.prog
                method = inf_obj.method
                ene = elstruct.reader.energy(prog, method, out_str)

                geo = elstruct.reader.opt_geometry(prog, out_str)

                save_path = tau_save_fs[-1].path(locs)
                print(" - Saving...")
                print(" - Save path: {}".format(save_path))

                tau_save_fs[-1].create(locs)
                tau_save_fs[-1].file.geometry_info.write(inf_obj, locs)
                tau_save_fs[-1].file.geometry_input.write(inp_str, locs)
                tau_save_fs[-1].file.energy.write(ene, locs)
                tau_save_fs[-1].file.geometry.write(geo, locs)

                # Saving the energy to a SP filesystem
                print(" - Saving energy of unique geometry...")
                sp_save_fs = autofile.fs.single_point(save_path)
                sp_save_fs[-1].create(mod_thy_info[1:4])
                sp_save_fs[-1].file.input.write(inp_str, mod_thy_info[1:4])
                sp_save_fs[-1].file.info.write(inf_obj, mod_thy_info[1:4])
                sp_save_fs[-1].file.energy.write(ene, mod_thy_info[1:4])

                saved_geos.append(geo)

        # update the tau trajectory file
        filesys.mincnf.traj_sort(tau_save_fs)


def assess_pf_convergence(tau_save_fs, ref_ene,
                          temps=(300., 500., 750., 1000., 1500.)):
    """ Determine how much the partition function has converged
    """

    # Calculate sigma values at various temperatures for the PF
    for temp in temps:
        sumq = 0.
        sum2 = 0.
        idx = 0
        print('integral convergence for T = ', temp)
        for locs in tau_save_fs[-1].existing():
            idx += 1
            ene = tau_save_fs[-1].file.energy.read(locs)
            ene = (ene - ref_ene) * phycon.EH2KCAL
            tmp = numpy.exp(-ene*349.7/(0.695*temp))
            sumq = sumq + tmp
            sum2 = sum2 + tmp**2
            sigma = numpy.sqrt(
                (abs(sum2/float(idx)-(sumq/float(idx))**2))/float(idx))
            print(sumq/float(idx), sigma, 100.*sigma*float(idx)/sumq, idx)


# eps[whatever], sig[ang] params
LJ_DCT = {
    ('H', 'H'): [0.25, 1.0],
    ('H', 'C'): [0.25, 1.0],
    ('H', 'O'): [0.25, 1.0],
    ('C', 'C'): [0.25, 1.0],
    ('C', 'O'): [0.25, 1.0],
    ('O', 'O'): [0.25, 1.0],
}

# A, B, C params E[kcal] R[Ang]
EXP6_DCT = {
    ('H', 'H'): [2.442e3, 3.74, 48.8],
    ('H', 'C'): [6.45e3, 3.67, 116.0],
    ('H', 'O'): [6.45e3, 3.67, 116.0],
    ('C', 'C'): [7.69e4, 3.6, 460.0],
    ('C', 'O'): [7.69e4, 3.6, 460.0],
    ('O', 'O'): [7.69e4, 3.6, 460.0]
}


def _low_repulsion_struct(zma_ref, zma_samp, thresh=20.0):
    """ Check if the coloumb sum
    """

    # # Convert to geoms
    geo_ref = automol.zmatrix.geometry(zma_ref)
    geo_samp = automol.zmatrix.geometry(zma_samp)

    # Calculate the pairwise potentials
    pot_mat = _pairwise_potential_matrix(geo_ref)
    pot_mat_samp = _pairwise_potential_matrix(geo_samp)

    # Generate the pairs for the potentials
    pairs = _generate_pairs(geo_ref)

    # Calculate sum of potentials
    sum_ref, sum_samp = 0.0, 0.0
    for (idx1, idx2) in pairs:
        sum_ref += pot_mat[idx1, idx2]
        sum_samp += pot_mat_samp[idx1, idx2]

    print('long_range_pots {:.2f} {:.2f} {:.2f}'.format(sum_ref, sum_samp, sum_samp-sum_ref))

    # # Check if the potentials are within threshold
    low_repulsion = bool((sum_samp - sum_ref) <= thresh)
    # low_repulsion = True

    return low_repulsion


def _generate_pairs(geo):
    """ Generate a list of pairs to calculate potentials
    """

    # Grab the indices of the heavy atoms for the zmas
    # heavy_idxs = automol.geom.atom_indices(geo, 'H', match=False)

    # Get the h atom idxs as list of list for each heavy atom
    # gra = automol.geom.graph(geo)
    # neigh_dct = automol.graph.atom_neighbor_keys(gra)
    # h_idxs = ()
    # for idx in heavy_idxs:
    #     neighs = neigh_dct[idx]
    #     h_idxs += (tuple(x for x in neighs if x not in heavy_idxs),)

    # Get heavy atom pairs
    # heavy_pairs = tuple(itertools.combinations(heavy_idxs, 2))
    # h_pairs = tuple()
    # for comb in itertools.combinations(h_idxs, 2):
    #     h_pairs += tuple(itertools.product(*comb))

    # Loop over neighbor dict to build pairs
    # print('geo', geo)
    # print('heavy idxs', heavy_idxs)
    # print('heacy pairs', heavy_pairs)
    # print('h idxs', h_idxs)
    # print('h pairs', h_pairs)
    # print(neigh_dct)
    #

    pairs = tuple()
    for i in range(len(geo)):
        for j in range(len(geo)):
            if i != j:
                pairs += ((i, j),)

    return pairs


def _pairwise_potential_matrix(geo):
    """ Generate a matrix of pairwise potentials
    """

    # Initialize matrix
    natoms = len(geo)
    pot_mat = numpy.zeros((natoms, natoms))

    for i in range(natoms):
        for j in range(natoms):
            pot_mat[i, j] = _pairwise_potentials(geo, (i, j))

    return pot_mat


def _pairwise_potentials(geo, idx_pair, potential='exp6'):
    """ Calculate the sum of the pairwise potential for a
        given set of atom pairs
    """

    # Get the indexes and symbols
    idx1, idx2 = idx_pair
    if idx1 != idx2:

        # Get the symbols of the atoms
        symbols = automol.geom.symbols(geo)
        symb1, symb2 = symbols[idx1], symbols[idx2]

        # Calculate interatomic distance
        rdist = automol.geom.distance(geo, idx1, idx2) * phycon.BOHR2ANG

        # Calculate the interaction potential value
        if potential == 'exp6':
            pot_val = _pairwise_exp6_potential(rdist, symb1, symb2)
        elif potential == 'lj_12_6':
            pot_val = _pairwise_lj_potential(rdist, symb1, symb2)
        else:
            pot_val = None

    else:
        pot_val = 1e10

    return pot_val


def _pairwise_exp6_potential(rdist, symb1, symb2):
    """ Calcualte pot
    """

    exp6_params = EXP6_DCT.get((symb1, symb2), None)
    if exp6_params is None:
        exp6_params = EXP6_DCT.get((symb2, symb1), None)

    pot_val = _exp6_potential(rdist, *exp6_params)

    return pot_val


def _exp6_potential(rdist, apar, bpar, cpar):
    """ Calculate modified Buckhingham potential
    """
    return apar * numpy.exp(-1.0*bpar*rdist) - (cpar / rdist**6)


def _pairwise_lj_potential(rdist, symb1, symb2):
    """ Calcualte pot
    """

    ljparams = LJ_DCT.get((symb1, symb2), None)
    if ljparams is None:
        ljparams = LJ_DCT.get((symb2, symb1), None)

    pot_val = _lj_potential(rdist, *ljparams)

    return pot_val


def _lj_potential(rdist, eps, sig):
    """ Calculate Lennard-Jones Potential
    """
    return (4.0 * eps) * ((sig / rdist)**12 - (sig / rdist)**6)
