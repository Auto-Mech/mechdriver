""" Generate the information necessary to product the vrctst input files
"""

import automol
from automol.graph import atom_sorted_neighbor_atom_keys


def intramolecular_constraint_dct(inf_sep_zma, rct_zmas):
    """ Set the additional constraints for the constrained MEP,
        constrains all of the intramolecular constraint_dct

        maybe just separate inf sep zma anf get the frag atoms
        sort of assumes the zma is frag1_zma + frag2_zma
    """

    frag1_natom = automol.zmat.count(rct_zmas[0])
    frag2_natom = automol.zmat.count(rct_zmas[1])

    # Build pairs for intermolecular coords to optimize:
    #   (zma_atom_idx, coord_idx in row) (uses 0-indexing)
    # frag1_natom, 0 is the scan coord already accounted for
    no_frz_idxs = []
    no_frz_idxs.append([frag1_natom, 0])
    no_frz_idxs.append([frag1_natom, 1])
    no_frz_idxs.append([frag1_natom, 2])
    if frag2_natom == 2:
        no_frz_idxs.append([frag1_natom+1, 1])
        no_frz_idxs.append([frag1_natom+1, 2])
    elif frag2_natom > 2:
        no_frz_idxs.append([frag1_natom+1, 1])
        no_frz_idxs.append([frag1_natom+1, 2])
        no_frz_idxs.append([frag1_natom+2, 1])

    # Now grab the coordinates NOT in the opt coord idxs
    alt_froz_coords = []
    name_matrix = automol.zmat.name_matrix(inf_sep_zma)
    for row_idx, row in enumerate(name_matrix):
        for coord_idx, coord in enumerate(row):
            if [row_idx, coord_idx] not in no_frz_idxs:
                if coord is not None:
                    alt_froz_coords.append(coord)

    # Now build the constraint dictionary
    zma_vals = automol.zmat.value_dictionary(inf_sep_zma)
    constraint_dct = dict(zip(
        alt_froz_coords, (zma_vals[name] for name in alt_froz_coords)))

    return constraint_dct


def fragment_geometries(ts_zma, rct_zmas, bnd_keys):
    """ Generates geometries of each reacting species where a dummy atom has
        been added at the location of an atom that will undergo bond formation
        with the reacting fragment.

        Notationally, for (R1-B1)---(B2--R2)

        where R-B1 and R2-B2 represent two reacting radicals where we
        have highlighted the Rn-Bn bond including the atom Bn which is involved
        in bond formation (B1-b2)'of the two radicals.
        We hope to generate geometriex

        Function takes a ZMA along the reaction MEP where bond formation is

        For example: For CH3, a dummy atom would be placed perpendicular
        to plane since any bond formed by CH3 will involve an atom located
        in this direction. For radicals, this often mimics the location of
        the reacting radical orbital.

        The location of this dummy atom is determined by fragmenting a reactive
        TS complex containing both fragments and severing it using the keys of
        the forming/breaking bond.

        Right now, it is assumed the atom ordering of each fragment in the
        ts_zma matches that of the reactants provided separately in rct_zmas.

        Also assumes ts_zma = frag1_zma + frag2_zma

        :param ts_zma: TS Z-Matrix containing reactive fragment
            at intermolecular distance where bond breaking/forming occuring
        :param rct_zmas: Z-Matrix for each reacting species alone, without
            influence of other fragment (species at infinite separation)
        :param bond-keys: atom indices for forming/breaking bond
    """

    def _get_chain_idx(geo, idx, excl_idxs=()):
        """ get neighbor to atom specified by idx. Take nonHs first
        """
        gra = automol.geom.graph(geo)
        chain_idxs = atom_sorted_neighbor_atom_keys(
            gra, idx, excl_keys=excl_idxs, symbs_last=('H',))
        return chain_idxs[0]

    # Group the forming bond where higher value is first (MAX, MIN)
    bnd_keys = sorted(list(bnd_keys), reverse=True)

    # Get the geometry and of the point on the MEP
    mep_total_geo = automol.zmat.geometry(ts_zma)
    mep_fgeos = [mep_total_geo[:bnd_keys[0]], mep_total_geo[bnd_keys[0]:]]

    # Get isolated fragments opt'd at infinite separation, aligned to MEP geoms
    iso_fgeos = [automol.zmat.geometry(zma) for zma in rct_zmas]
    (iso1_symbs, iso2_symbs) = (automol.geom.symbols(geo) for geo in iso_fgeos)
    (mep1_symbs, mep2_symbs) = (automol.geom.symbols(geo) for geo in mep_fgeos)
    if iso1_symbs != mep1_symbs or iso2_symbs != mep2_symbs:
        iso_fgeos[0], iso_fgeos[1] = iso_fgeos[1], iso_fgeos[0]
    iso_fgeos = [automol.geom.align(geo, mep_fgeos[i])
                 for i, geo in enumerate(iso_fgeos)]

    # Now build the list of the indices for the following loops
    # Used to define the X via set of internal coordinates
    # x_idx: index for geom to place dummy X atom
    # a1_idx: index corresponding to "bonding" atom in geometry
    # a2_idx and a3_idx are just some atom down a chain, preferably nonHs
    # Because of the align procedure we can assume the idxs in MEP and iso same
    # possible this will only work for H (with regard to 2nd set of frag idxs)
    frag1_natoms = len(mep_fgeos[0])
    x_idx = frag1_natoms
    a1_idx = bnd_keys[1]
    a2_idx = _get_chain_idx(iso_fgeos[0], a1_idx, excl_idxs=())
    if frag1_natoms > 3:
        a3_idx = _get_chain_idx(iso_fgeos[0], a2_idx, excl_idxs=(a1_idx,))
    elif frag1_natoms == 3:
        # has to be 0,1,2; for example: if a1=1, a2=0 then a3=2
        a3_idx = next(iter({0, 1, 2} - {a1_idx, a2_idx}))
    else:
        a3_idx = None
    frag_idxs = (
        (x_idx, a1_idx, a2_idx, a3_idx),
        (0, 1, 2, 3)
    )

    # Add dummy to each isolated frag geom to simulate atom its bound to in TS
    # Need to generate coordinates of dummy in isolated frag geom system by
    # (1) Calculate values of internal coords that define X
    #     relative to 4 atoms in the MEP geom
    # (2) Use values from (1) and frag xyzs to calculate xyz of X in frag geom
    iso_fgeos_wdummy = ()
    mol_data = zip(mep_fgeos, iso_fgeos, bnd_keys, frag_idxs)
    for i, (mep_fgeo, iso_fgeo, bnd_key, fidxs) in enumerate(mol_data):

        if not automol.geom.is_atom(mep_fgeo):

            # Build MEPFragGeom+X coordinates using MEP geometry
            # x_coord defined by frm_keys
            x_coord = mep_total_geo[bnd_key][1]
            dummy_row = ('X', x_coord)
            if i == 0:
                mep_geox = mep_fgeo + (dummy_row,)
            else:
                mep_geox = (dummy_row,) + mep_fgeo

            # Calculate coords to define X position in IsoFragGeom structure
            xdistance = automol.geom.distance(mep_geox, *fidxs[0:2])
            xangle = automol.geom.central_angle(mep_geox, *fidxs[0:3])
            if len(mep_fgeo) > 2:
                xdihedral = automol.geom.dihedral_angle(mep_geox, *fidxs[0:4])
            else:
                xdihedral = 0.0

            # Set the coords for the IsoFragStructure
            xyz1 = iso_fgeo[fidxs[1]][1]
            xyz2 = iso_fgeo[fidxs[2]][1]
            xyz3 = iso_fgeo[fidxs[3]][1] if len(mep_fgeo) > 2 else 0.0

            # Calculate the X Position for the IsoFrag structure
            xyzp = automol.util.vector.from_internals(
                dist=xdistance, xyz1=xyz1, ang=xangle, xyz2=xyz2,
                dih=xdihedral, xyz3=xyz3)

            # Generate the IsoFragGeom+X coordinates for the structure.inp file
            if i == 0:
                iso_geo_wdummy = iso_fgeo + (('X', xyzp),)
            else:
                iso_geo_wdummy = (('X', xyzp),) + iso_fgeo

            # Append to final geoms
            iso_fgeos_wdummy += (iso_geo_wdummy,)

        else:
            # If atom, set IsoFragGeom+X coords equal to mep_geo
            iso_fgeos_wdummy += (mep_fgeo,)

    return mep_total_geo, iso_fgeos_wdummy, (frag_idxs[0][1], frag_idxs[1][1])


# def assess_face_symmetries(divsur_out_string):
def assess_face_symmetries(fgeo1, fgeo2):
    """ check the symmetry of the faces for each fragment

        :param fgeo1: fragment geometry 1
        :type fgeo1: molecular geometry 1
    """

    # Read fragment geoms from divsur.out with coordinates in the divsur frame
    # fgeo1, fgeo2 = varecof_io.reader.divsur.frag_geoms_divsur_frame(
    #     divsur_out_string)
    # fgeos = [automol.geom.from_string(fgeo1),
    #          automol.geom.from_string(fgeo2)]

    # Check facial symmetry if fragments are molecules
    symms = [False, False]
    # for i, fgeo in enumerate(fgeos):
    for i, fgeo in enumerate((fgeo1, fgeo2)):
        if not automol.geom.is_atom(fgeo):
            # Reflect the dummy atom (pivot position) about the xy plane
            if i == 0:
                dummy_idx = len(fgeo) - 1
            else:
                dummy_idx = 0
            fgeo_reflect = automol.geom.reflect_coordinates(
                fgeo, [dummy_idx], ['x', 'y'])
            # Compute Coloumb spectrum for each geom to its reflected version
            symms[i] = automol.geom.almost_equal_coulomb_spectrum(
                fgeo, fgeo_reflect, rtol=5e-2)

    # Set the face and face_sym keywords based on the above tests
    # [symm1, symm2] = symms
    # if symm1 and symm2:
    #     faces = [0, 1]
    #     face_symm = 4
    # elif symm1 and not symm2:
    #     faces = [0, 1]
    #     face_symm = 2
    # elif not symm1 and symm2:
    #     faces = [0, 1]
    #     face_symm = 2
    # elif not symm1 and not symm2:
    #     faces = [0]
    #     face_symm = 1

    faces = [0]
    face_symm = 1
    return faces, face_symm


# FUNCTIONS TO SET UP THE DIVIDING SURFACE FRAMES
def build_pivot_frames(frag_geos_wdummy, frag_a1_idxs):
    """ Use geometries to get pivot info only set up for 1 or 2 pivot points

        change to just pass the bonding atom idx (bonding atom, the one that
        bonds to dummy)
        use the assumption about the a1 idx used in the fragment geoms?
    """

    frames, npivots, = (), ()
    for i, (geo, a1_idx) in enumerate(zip(frag_geos_wdummy, frag_a1_idxs)):

        geom = automol.geom.without_dummy_atoms(geo)

        # Single pivot point centered on atom
        if automol.geom.is_atom(geom):
            # Single pivot point centered on atom
            npivot = 1
            frame = (0, 0, 0, 0)
        elif automol.geom.is_diatomic(geom):
            # For linear species we place the pivot point on radical
            # with no displacment, so no need to find coordinates
            npivot = 1
            frame = (0, 0, 0, 0)
        else:
            # else we build an xy frame to easily place pivot point
            # explain the frame???
            npivot = 2

            # For each fragment, get indices for a
            # chain (up to three atoms, that terminates at the dummy atom)
            gra = automol.geom.graph(geom)
            gra_neighbor_dct = automol.graph.atoms_neighbor_atom_keys(gra)
            bond_neighbors = gra_neighbor_dct[a1_idx]

            # Find idx in each fragment geom that corresponds to the bond index
            for j, idx in enumerate(bond_neighbors):
                if geom[idx][0] != 'H':
                    bond_neighbor_idx = idx
                    break
                if geom[idx][0] == 'H' and j == (len(bond_neighbors) - 1):
                    bond_neighbor_idx = idx

            # Set up the frame indices for the divsur file
            if i == 0:
                pivot_idx = len(geom)
                frame = (a1_idx, pivot_idx, bond_neighbor_idx, a1_idx)
            else:
                pivot_idx = 0
                a1_idx += 1
                bond_neighbor_idx += 1
                frame = (a1_idx, pivot_idx, bond_neighbor_idx, a1_idx)
            frame = tuple(val+1 for val in frame)

        # Append to lists
        frames += (frame,)
        npivots += (npivot,)

    return frames, npivots


def calc_pivot_angles(frag_geos_wdummy, frames):
    """ Calculate the angle for three atoms which define the frame
        of the dividing surface.

        :param frag_geos_wdummy: geometries of fragments at infinite
            separation, including dummy atoms showing location of
            reactive atom from other fragment
        :type frag_geos_wdummy: tuple(automol.geom object)
        :param frames: 4-idx set for the atoms in each fragment that
            define the orientational axes of the dividing surface frame.
        :type frames: tuple(tuple(int))
    """

    angles = tuple()
    for geo_wdummy, frame in zip(frag_geos_wdummy, frames):
        geo_ndum = automol.geom.without_dummy_atoms(geo_wdummy)
        if automol.geom.is_atom(geo_ndum) or automol.geom.is_linear(geo_ndum):
            angle = None
        else:
            frame = [val-1 for val in frame]
            angle = automol.geom.central_angle(
                geo_wdummy, frame[2], frame[0], frame[1])

        angles += (angle,)

    return angles


def calc_pivot_xyzs(total_geo, frag_geos, bnd_keys):
    """ Determine the fragment xyz-coordinates that will be used
        to set the position of the pivot points.

        We only need to determine the coord for diatomics, which
        are set to the reactive radical site. Atoms and polyatomics
        will be set using a special frame where the origin is
        sufficient.

        :param total_geo: geom with both reactive fragments
            where bond-formation is occuring
        :type total_geo: automol.geom object
        :param frag_geos: geom of each fragment at infinite separation
        :type frag_geos: tuple(automol.geom object)
        :param bnd_keys: keys of reactive atoms in total_geo
        :type bnd_keys:
    """

    bnd_keys = sorted(list(bnd_keys))

    xyzs = tuple()
    for rxn_idx, geo in zip(bnd_keys, frag_geos):
        if automol.geom.is_diatomic(geo):
            xyz = total_geo[rxn_idx][1]
        else:
            xyz = (0.0, 0.0, 0.0)

        xyzs += (xyz,)

    return xyzs
