""" Generate the information necessary to product the vrctst input files
"""

import os
import subprocess
import automol
import varecof_io


def calc_vrctst_rates():
    """ From scripts.es
    """

    # Set paths and build dirs for VRC-TST calculation is run
    vrc_path = os.path.join(os.getcwd(), 'vrc')
    scr_path = os.path.join(vrc_path, 'scratch')
    os.makedirs(vrc_path, exist_ok=True)
    os.makedirs(scr_path, exist_ok=True)
    print('vrc_path test:', vrc_path)

    # Correction potential
    corr_pot = False
    if corr_pot:
        # Read the values for the correction potential from filesystem
        potentials, pot_labels = moldr.vrctst.read_corrections()
        
        # Build correction potential .so file used by VaReCoF
        moldr.vrctst. build_correction_potential(
            mep_distances, potentials,
            bnd_frm_idxs, fortran_compiler, vrc_path,
            dist_restrict_idxs=(),
            pot_labels=pot_labels,
            pot_file_names=(),
            spc_name='mol')

    # Write the electronic structure template file
    memory = 4.0
    basis = 'cc-pvdz'
    method = '{rs2c, shift=0.25}'

    num_act_elc = high_mul
    num_act_orb = num_act_elc
    ts_formula = automol.geom.formula(automol.zmatrix.geometry(ts_zma))
    _, wfn_str = moldr.ts.cas_options_2(
        ts_info, ts_formula, num_act_elc, num_act_orb, high_mul)
    tml_inp_str = varecof_io.writer.input_file.tml(
        memory, basis, wfn_str, method, inf_sep_ene)

    # Write machines file to set compute nodes
    machines = ['b450:8', 'b451:8', 'b452:8', 'b453:8']
    with open(os.path.join(vrc_path, 'machines'), 'w') as machine_file:
        for machine in machines:
            machine_file.writelines(machine + '\n')

    # Write the remaining input file strings
    input_strs = moldr.vrctst.input_prep(ts_zma, ts_dct['rct_zmas'], dist_name, vrc_path)
    [struct_inp_str, lr_divsur_inp_str, tst_inp_str,
     els_inp_str, mc_flux_inp_str, conv_inp_str] = input_strs
    
    with open(os.path.join(vrc_path, 'structure.inp'), 'w') as inp_file:
        inp_file.write(struct_inp_str)
    with open(os.path.join(vrc_path, 'lr_divsur.inp'), 'w') as inp_file:
        inp_file.write(lr_divsur_inp_str)
    with open(os.path.join(vrc_path, 'tst.inp'), 'w') as inp_file:
        inp_file.write(tst_inp_str)
    with open(os.path.join(vrc_path, 'molpro.inp'), 'w') as inp_file:
        inp_file.write(els_inp_str)
    with open(os.path.join(vrc_path, 'mc_flux.inp'), 'w') as inp_file:
        inp_file.write(mc_flux_inp_str)
    with open(os.path.join(vrc_path, 'convert.inp'), 'w') as inp_file:
        inp_file.write(conv_inp_str)
    with open(os.path.join(vrc_path, 'mol.tml'), 'w') as tml_file:
        tml_file.write(tml_inp_str)

    print('\n\nEXITING VRCTST AFTER WRITING INPUT')




def input_prep(ts_zma, rct_zmas, dist_name, vrc_path):
    """ prepare all the input files for a vrc-tst calculation
    """

    # Set info on the form indices
    bnd_frm_idxs = automol.zmatrix.bond_idxs(ts_zma, dist_name)
    min_idx, max_idx = min(bnd_frm_idxs), max(bnd_frm_idxs)

    # Build geometries needed for the varecof run
    total_geom, frag_geoms, frag_geoms_wdummy = fragment_geometries(
        ts_zma, rct_zmas, min_idx, max_idx)

    # Set information for the pivot points needed in divsur.inp
    frames, npivots = build_pivot_frames(
        min_idx, max_idx, total_geom, frag_geoms, frag_geoms_wdummy)
    pivot_angles = calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames)
    pivot_xyzs = calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms)

    # Write the long- and short-range divsur input files
    r1dists_lr = [8., 6., 5., 4.5, 4.]
    lr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_lr, 1, 1, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])

    # Write the short-range divsur files
    r1dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    r2dists_sr = [4., 3.8, 3.6, 3.4, 3.2, 3., 2.8, 2.6, 2.4, 2.2]
    d1dists = [0.01, 0.5, 1.]
    d2dists = [0.01, 0.5, 1.]
    t1angs = [pivot_angles[0]] if pivot_angles[0] is not None else []
    t2angs = [pivot_angles[1]] if pivot_angles[0] is not None else []
    if automol.geom.is_atom(frag_geoms[0]):
        d1dists = []
        t1angs = []
    if automol.geom.is_atom(frag_geoms[1]):
        d2dists = []
        t2angs = []
    if automol.geom.is_linear(frag_geoms[0]):
        d1dists = [0.]
        t1angs = []
    if automol.geom.is_linear(frag_geoms[1]):
        d2dists = [0.]
        t2angs = []
    conditions = {}
    # conditions = {'delta_r': 0}
    sr_divsur_inp_str = varecof_io.writer.input_file.divsur(
        r1dists_sr, npivots[0], npivots[1], pivot_xyzs[0], pivot_xyzs[1],
        frame1=frames[0],
        frame2=frames[1],
        d1dists=d1dists,
        d2dists=d2dists,
        t1angs=t1angs,
        t2angs=t2angs,
        r2dists=r2dists_sr,
        **conditions)

    # Write the structure input files
    struct_inp_str = varecof_io.writer.input_file.structure(
        frag_geoms_wdummy[0], frag_geoms_wdummy[1])
    with open(os.path.join(vrc_path, 'structure.inp'), 'w') as inp_file:
        inp_file.write(struct_inp_str)

    # Write the divsur input file with determined frames
    with open(os.path.join(vrc_path, 'divsur.inp'), 'w') as inp_file:
        inp_file.write(sr_divsur_inp_str)

    # Obtain the divsur.out file with divsur-frame fragment geoms
    divsur_out_str = build_divsur_out_file(vrc_path, os.getcwd())

    # Write the tst.inp file
    nsamp_max = 2000
    nsamp_min = 50
    flux_err = 10
    pes_size = 2
    faces, faces_symm = assess_face_symmetries(divsur_out_str)
    tst_inp_str = varecof_io.writer.input_file.tst(
        nsamp_max, nsamp_min, flux_err, pes_size,
        faces=faces, faces_symm=faces_symm)

    # Write the potential energy surface input string
    exe_path = '/blues/gpfs/home/sjklipp/bin/molpro'
    base_name = 'mol'
    els_inp_str = varecof_io.writer.input_file.elec_struct(
        exe_path, base_name)

    # Write the mc_flux.inp input string
    mc_flux_inp_str = varecof_io.writer.input_file.mc_flux()

    # Write the convert.inp input string
    conv_inp_str = varecof_io.writer.input_file.convert()

    input_strs = [
        struct_inp_str, lr_divsur_inp_str, tst_inp_str,
        els_inp_str, mc_flux_inp_str, conv_inp_str]

    return input_strs


def read_potentials(
    rxn_save_fs, 
    mep_opt_thy_level,
    high_ene_thy_levels,
    mep_r_constrained_coords,
    mep_intramol_constrained_coords,
    r_grid_vals):
    """ Read values form the filesystem to get the values to
        correct ht MEP
    """

    # Hard part is the filesystem storage
    # FS
    # SCANS
    # Rn
    # 2.1, 2.4, 2.7, 3.0, 3.3, 3.6

    # larger scans:
    # scan over Rn
    # one dir has JUST Rn constrained
    # other has Rn constrained and others

    # accept as args, filesystems and/or thy levels for:
    # energies along MEP with full opts (except Rn rxn coord)
    # energies along MEP with opts of intermolecular coords
    # higher level energies along MEP 
   
    # calc relaxation correction
    # calc energy corrections for each higher level energy
    # calc combined relaxation-energy corrections
    
    # assume you have fs vars passed in 
    rxn_save_path = 'save_fs_arg'

    # Filesystem for all the optimizations (full, constrained) along MEP
    mep_opt_thy_save_fs = autofile.fs.theory(rxn_save_path)
    mep_opt_thy_save_fs.leaf.create(mep_opt_thy_level[1:4])
    mep_opt_scn_save_fs = autofile.fs.scan(mep_opt_thy_save_fs)        
    mep_opt_scn_save_fs_path = mep_opt_scn_save_fs.leaf.path() 

    # Read the energies from the full and constrained opts along MEP
    # for idx, 
    # locs = [[], []] 
    #     # get energy
    #     if not scn_save_fs.leaf.file.energy.exists(locs):
    #         continue
    #     else:
    #         ene = scn_save_fs.leaf.file.energy.read(locs)

    return potentials, potential_labels


def build_correction_potential(mep_distances, potentials,
                               bnd_frm_idxs, fortran_compiler, vrc_path,
                               dist_restrict_idxs=(),
                               pot_labels=(),
                               pot_file_names=(),
                               spc_name='mol'):
    """  use the MEP potentials to compile the correction potential .so file
    """

    # Read the potentials needed for corrections

    # Build string Fortan src file containing correction potentials
    species_corr_str = varecof_io.writer.corr_potentials.species(
        mep_distances,
        potentials,
        bnd_frm_idxs,
        dist_restrict_idxs=dist_restrict_idxs,
        pot_labels=pot_labels,
        species_name=spc_name)

    # Build string dummy corr file where no correction used
    dummy_corr_str = varecof_io.writer.corr_potentials.dummy()
    
    # Build string for auxiliary file needed for spline fitting
    pot_aux_str = varecof_io.writer.corr_potentials.auxiliary()

    # Build string for makefile to compile corr pot file into .so file
    makefile_str = varecof_io.writer.corr_potentials.makefile(
        fortran_compiler,
        pot_file_names=pot_file_names)

    # Write all of the files needed to build the correction potential
    with open(os.path.join(vrc_path, spc_name+'_corr.f'), 'w') as corr_file:
        corr_file.write(species_corr_str)
    with open(os.path.join(vrc_path, 'dummy_corr.f'), 'w') as corr_file:
        corr_file.write(dummy_corr_str)
    with open(os.path.join(vrc_path, 'pot_aux.f'), 'w') as corr_file:
        corr_file.write(pot_aux_str)
    with open(os.path.join(vrc_path, 'makefile'), 'w') as corr_file:
        corr_file.write(makefile_str)

    # Compile the correction potential
    varecof_io.writer.corr_potentials.compile_corr_pot(vrc_path)


def fragment_geometries(ts_zma, rct_zmas, min_idx, max_idx):
    """ Generate the fragment geometries from the ts Z-matrix and the
        indices involved in the forming bond
    """

    # Get geometries of fragments from the ts_zma from the MEP
    mep_total_geo = automol.zmatrix.geometry(ts_zma)
    mep_fgeos = [mep_total_geo[:max_idx], mep_total_geo[max_idx:]]

    # Get geometries of isolated fragments at infinite sepearation
    iso_fgeos = [automol.zmatrix.geometry(zma) for zma in rct_zmas]

    # displace geoms to avoid having 0.0s
    # iso_fgeos = [automol.geom.translated(geo, (10.0, 10.0, 10.0)) for geo in iso_fgeos]

    # Hack to get my rct_zmas in right order for one example. Need better fix
    iso_fgeos[0], iso_fgeos[1] = iso_fgeos[1], iso_fgeos[0]

    # Get the geometries for the structure.inp file
    iso_fgeos_wdummy = []
    mol_data = zip(mep_fgeos, iso_fgeos, (max_idx, min_idx))
    for i, (mep_fgeo, iso_fgeo, idx) in enumerate(mol_data):

        if not automol.geom.is_atom(mep_fgeo):

            # Build MEPFragGeom+X coordinates using MEP geometry
            # x_idx: index for geom to place dummy X atom
            # a1_idx: index corresponding to "bonding" atom in geometry
            x_coord = mep_total_geo[idx][1]
            dummy_row = ('X', x_coord)
            if i == 0:
                mep_geo_wdummy = mep_fgeo + (dummy_row,)
                x_idx = len(mep_geo_wdummy) - 1
                a1_idx = 0
            else:
                mep_geo_wdummy = (dummy_row,) + mep_fgeo
                x_idx = 0
                a1_idx = 1

            # Set a2_idx to a1_idx + 1 as I don't think there are any restrictions
            a2_idx = a1_idx + 1

            # Set a3_idx. 
            # Need to ensure idx does NOT correspond to atom where x = 0.0
            # The internal xyzp routine dies in this case
            for idx in range(a2_idx+1, len(iso_fgeo)):
                if not iso_fgeo[idx][1][0] == 0.0:
                    a3_idx = idx
                    break

            # Calculate coords to define X position in IsoFragGeom structure
            print('iso_fgeo')
            print(iso_fgeo)
            print('mep_geo')
            print(automol.geom.string(mep_geo_wdummy))
            print('idxs')
            print(x_idx)
            print(a1_idx)
            print(a2_idx)
            print(a3_idx)
            xyz1 = iso_fgeo[a1_idx][1]
            xyz2 = iso_fgeo[a2_idx][1]
            xdistance = automol.geom.distance(
                mep_geo_wdummy, x_idx, a1_idx)
            xangle = automol.geom.central_angle(
                mep_geo_wdummy, x_idx, a1_idx, a2_idx)
            if len(mep_fgeo) > 2:
                xyz3 = iso_fgeo[a3_idx][1]
                xdihedral = automol.geom.dihedral_angle(
                    mep_geo_wdummy, x_idx, a1_idx, a2_idx, a3_idx)
            else:
                xyz3 = 0.0
                xdihedral = 0.0

            # Calculate the X Position for the IsoFrag structure
            print('\n\ninternals info')
            print(xyz1)
            print(xyz2)
            print(xyz3)
            print(xdistance)
            print(xangle)
            print(xdihedral)
            xyzp = automol.geom.find_xyzp_using_internals(
                xyz1, xyz2, xyz3, xdistance, xangle, xdihedral)

            # Generate the IsoFragGeom+X coordinates for the structure.inp file
            if i == 0:
                iso_geo_wdummy = iso_fgeo + (('X', xyzp),)
            else:
                iso_geo_wdummy = (('X', xyzp),) + iso_fgeo

            # Append to final geoms
            iso_fgeos_wdummy.append(iso_geo_wdummy)

        else:
            # If atom, set IsoFragGeom+X coords equal to mep_geo
            iso_fgeos_wdummy.append(mep_fgeo)

    return mep_total_geo, iso_fgeos, iso_fgeos_wdummy


def build_divsur_out_file(vrc_path, work_path):
    """ get the divsur.out string containing divsur-frame geoms
    """
    
    # Have to to path with divsur.inp to run script (maybe can fix)
    os.chdir(vrc_path)

    # Run the VaReCoF utility script to get the divsur.out file
    # Contains the fragment geometries in the divsur-defined coord sys
    varecof_io.writer.util.divsur_frame_geom_script()

    # Read fragment geoms from divsur.out with coordinates in the divsur frame
    with open(os.path.join(vrc_path, 'divsur.out'), 'r') as divsur_file:
        output_string = divsur_file.read()

    os.chdir(work_path)

    return output_string


def assess_face_symmetries(divsur_out_string):
    """ check the symmetry of the faces for each fragment
    """

    # Read fragment geoms from divsur.out with coordinates in the divsur frame
    fgeo1, fgeo2 = varecof_io.reader.divsur.frag_geoms_divsur_frame(
        divsur_out_string)
    fgeos = [
        automol.geom.from_string(fgeo1),
        automol.geom.from_string(fgeo2)
    ]

    # Check facial symmetry if fragments are molecules
    symms = []
    for i, fgeo in enumerate(fgeos):

        if not automol.geom.is_atom(fgeo):
            # Reflect the dummy atom (pivot position) about the xy plane
            if i == 0:
                dummy_idx = len(fgeo) - 1
            else:
                dummy_idx = 0
            fgeo_reflect = automol.geom.reflect_coordinates(
                fgeo, [dummy_idx], ['x', 'y'])
            # Compute Coloumb spectrum for each geom to its reflected version
            symm = automol.geom.almost_equal_coulomb_spectrum(
                fgeo, fgeo_reflect, rtol=5e-2)
        else:
            symm = False

        symms.append(symm)

    # Set the face and face_sym keywords based on the above tests
    [symm1, symm2] = symms
    if symm1 and symm2:
        faces = [0, 1]
        face_symm = 4
    elif symm1 and not symm2:
        faces = [0, 1]
        face_symm = 2
    elif not symm1 and symm2:
        faces = [0, 1]
        face_symm = 2
    elif not symm1 and not symm2:
        faces = [0]
        face_symm = 1

    return faces, face_symm


def build_pivot_frames(min_idx, max_idx,
                       total_geom, frag_geoms, frag_geoms_wdummy):
    """ use geometries to get pivot info
        only set up for 1 or 2 pivot points
    """

    frames, npivots, = [], []
    geom_data = zip([min_idx, max_idx], frag_geoms, frag_geoms_wdummy)
    for i, (rxn_idx, geom, geom_wdummy) in enumerate(geom_data):

        # Single pivot point centered on atom
        print('geom test in build_pivot:', geom)
        if automol.geom.is_atom(geom):
            npivot = 1
            frame = [0, 0, 0, 0]
        # For linear species we place the pivot point on radical
        # with no displacment, so no need to find coordinates
        elif automol.geom.is_linear(geom):
            npivot = 2
            frame = [0, 0, 0, 0]
        else:
            # else we build an xy frame to easily place pivot point
            npivot = 2

            # Find the idx in each fragment bonded to the atom at the pivot pt
            for j, coords in enumerate(geom):
                if coords == total_geom[rxn_idx]:
                    bond_idx = j
                    break

            # For each fragment, get indices for a
            # chain (up to three atoms, that terminates at the dummy atom)
            gra = automol.geom.graph(geom)
            gra_neighbor_dct = automol.graph.atom_neighbor_keys(gra)
            bond_neighbors = gra_neighbor_dct[bond_idx]

            # Find idx in each fragment geom that corresponds to the bond index
            for j, idx in enumerate(bond_neighbors):
                if geom[idx][0] != 'H':
                    bond_neighbor_idx = idx
                    break
                elif geom[idx][0] == 'H' and j == (len(bond_neighbors) - 1):
                    bond_neighbor_idx = idx

            # Set up the frame indices for the divsur file
            if i == 0:
                pivot_idx = len(geom)
                frame = [bond_idx, bond_neighbor_idx, pivot_idx, bond_idx]
            else:
                pivot_idx = 0
                bond_idx += 1
                bond_neighbor_idx += 1
                frame = [bond_idx, bond_neighbor_idx, pivot_idx, bond_idx]
            frame = [val+1 for val in frame]

        # Append to lists
        frames.append(frame)
        npivots.append(npivot)

    return frames, npivots


def calc_pivot_angles(frag_geoms, frag_geoms_wdummy, frames):
    """ get the angle for the three atoms definining the frame
    """
    angles = []
    for geom, geom_wdummy, frame in zip(frag_geoms, frag_geoms_wdummy, frames):
        print('pivot angles info')
        print(geom)
        print(geom_wdummy)
        print(frame)
        if automol.geom.is_atom(geom) or automol.geom.is_linear(geom):
            angle = None
        else:
            frame = [val-1 for val in frame]
            angle = automol.geom.central_angle(
                geom_wdummy, frame[2], frame[0], frame[1])

        angles.append(angle)

    return angles


def calc_pivot_xyzs(min_idx, max_idx, total_geom, frag_geoms):
    """ figure out where pivot point will be centered
        only linear speces need to have non-zero xyz, as they
        will not have a frame set up for them like atoms and
        polyatomics
    """
    xyzs = []
    for rxn_idx, geom in zip([min_idx, max_idx], frag_geoms):
        if automol.geom.is_linear(geom):
            xyz = total_geom[rxn_idx][1]
        else:
            xyz = [0.0, 0.0, 0.0]

        xyzs.append(xyz)

    return xyzs
