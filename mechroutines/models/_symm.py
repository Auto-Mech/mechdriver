""" Handle symmetry factor stuff
"""

import numpy

import automol
from autofile import fs
from mechlib.amech_io import printer as ioprinter

def symmetry_factor(pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors,
                    grxn=None, zma=None):
    """ Determines the the overall (internal and external) symmetry factor for
        a species or saddle point.

        Function will simply take a symmetry factor provided by the user by way
        of the spc_mod_dct_i, else it will calculate the symmetry factor using
        the requested procedure.

        For saddle points, the function ignores the possibility that two
        configurations differ only in their torsional values. As a result,
        the symmetry factor is a lower bound of the true value.

        :param pf_filesystems:
        :param grxn:
        :rtype: float
    """

    symm_factor = spc_dct_i.get('symm_factor')
    if symm_factor is not None:
        ioprinter.info_message(
            ' - Reading symmetry number input by user:', symm_factor)
    else:

        zrxn = spc_dct_i.get('zrxn', None)
        if zrxn is not None:
            grxn = automol.reac.relabel_for_geometry(zrxn)
        else:
            grxn = None

        sym_model = spc_mod_dct_i['symm']['mod']

        # Obtain geometry, energy, and symmetry filesystem

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':
            [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['symm']
            geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
            # Obtain the external symssetry number
            ext_symm = automol.geom.external_symmetry_factor(geo)

            # Set up the symmetry filesystem, read symmetrically similar geos
            # includes minimum geo
            sym_fs = fs.symmetry(cnf_path)
            symm_geos = [geo]
            symm_geos += [sym_fs[-1].file.geometry.read(locs)
                          for locs in sym_fs[-1].existing()]

            # Obtain the internal symmetry number and end group factors
            if rotors is not None:
                ioprinter.info_message(
                    ' - Determining internal sym number ',
                    'using sampling routine.')
                int_symm, endgrp = automol.symm.internal_symm_from_sampling(
                    symm_geos, rotors, grxn=grxn, zma=zma)
            else:
                ioprinter.info_message(' - No torsions, internal sym is 1.0')
                int_symm, endgrp = 1.0, 1.0

            # Obtain overall number, reduced as needed
            int_symm = automol.symm.reduce_internal_symm(
                geo, int_symm, ext_symm, endgrp)

        elif sym_model == 'HCO_model':
            if zma is not None:
                geo = automol.zmat.geometry(zma)
            else:
                [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['symm']
                geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
            ret = automol.symm.oxygenated_hydrocarbon_symm_num(geo, grxn)
            int_symm, ext_symm = ret

        else:
            [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['symm']
            geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)
            ioprinter.info_message(
                'No symmetry model requested, ',
                'setting internal sym factor to 1.0')
            ext_symm = automol.geom.external_symmetry_factor(geo)
            int_symm = 1.0

        if rotors is not None:
            rotor_symms = automol.rotor.symmetries(rotors, flat=True)
            int_symm = automol.symm.rotor_reduced_symm_factor(
                int_symm, rotor_symms)

            ext_symm *= _umbrella_factor(rotors, geo)

        symm_factor = ext_symm * int_symm

    return symm_factor


def _umbrella_factor(rotors, geo, grxn=None):
    """ check to see if this torsion has umbrella floppies
    """
    umb_fact = 1
    gra = automol.graph.dominant_resonance(automol.geom.graph(geo))
    rad_atms = automol.graph.radical_atom_keys(gra)
    # dih check
    adj_atms_dct = automol.graph.atoms_neighbor_atom_keys(gra)
    for rad_atm in rad_atms:
        adj_atms = adj_atms_dct[rad_atm]  
        if len(adj_atms) == 3:
            a, b, c = adj_atms
            dih_ang = automol.geom.dihedral_angle(geo, a, rad_atm, b, c)
            planarity = min([abs(dih_ang + x - numpy.pi) for x in [-2*numpy.pi, 0, 2*numpy.pi, 4*numpy.pi]])
            print(' umbrella mode dihedral', planarity)
        else:
            print('no umbrel dihedral')
    if rad_atms:
        dont_dbl = ()
        for rotor in rotors:
            for torsion in rotor:
                if any([ax_atm in dont_dbl for ax_atm in torsion.axis]):
                    continue
                if any([rad_atm in torsion.axis for rad_atm in rad_atms]):
                    if automol.pot.is_symmetric(torsion.pot):
                        dont_dbl += tuple([rad_atm for rad_atm in torsion.axis if rad_atm in rad_atms])
                        umb_fact *= 2
                        ioprinter.info_message(
                            'Umbrella mode identified for torsion about', torsion.axis)
    return umb_fact
