""" Handle symmetry factor stuff
"""

import numpy

import automol
from autofile import fs
from mechlib.amech_io import printer as ioprinter

def symmetry_factor(pf_filesystems, spc_mod_dct_i, spc_dct_i, rotors,
                    grxn=None, zma=None, racemic=True):
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
            zrxn = automol.reac.without_structures(zrxn)
            grxn = automol.reac.undo_zmatrix_conversion(zrxn)
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
                int_symm, endgrp = automol.symm.symmetry_factors_from_sampling(
                    symm_geos, rotors, grxn=grxn)
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
            ret = automol.symm.oxygenated_hydrocarbon_symm_num(geo, grxn, racemic=racemic)
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
            rotor_symms = automol.data.rotor.rotors_torsion_symmetries(rotors, flat=True)
            int_symm = automol.symm.rotor_reduced_symm_factor(
                int_symm, rotor_symms)

            ext_symm *= _umbrella_factor(rotors, geo)
        # if spc_dct_i['smiles'] in ['CC(OO)C[CH2]', '[CH2]CCCOO','CC(OO)CC[CH2]','[CH2]CCCCOO','CC(C[CH2])COO','OOCC([CH2])(C)C','OOCC([CH2])C','OOCCC(C)(C)[CH2]','OOC(C([CH2])(C)C)C','OOCC(C([CH2])C)C','CCC(COO)([CH2])C','OOC(C([CH2])C)C(C)C','OOC(CC(C)C)([CH2])C','OOC(CC([CH2])(C)C)(C)C','OOC(C(C)(C)C)C(C)[CH2]','OOCC(CC(C)C)([CH2])C','OOC(CC(C)(C)C)([CH2])C','OOC(C([CH2])(C)C)C(C)C']:
        #     print('divide ext sym by two')
        #     ext_symm /= 2
        symm_factor = ext_symm * int_symm

    return symm_factor


def _umbrella_factor(rotors, geo, grxn=None):
    """ check to see if this torsion has umbrella floppies
    """
    gra = automol.graph.kekule(automol.geom.graph(geo))
    rad_atms = automol.graph.radical_atom_keys(gra)
    # dih check
    adj_atms_dct = automol.graph.atoms_neighbor_atom_keys(gra)
    planarity = 100.
    umb_inv = False
    for rad_atm in rad_atms:
        adj_atms = adj_atms_dct[rad_atm]
        if len(adj_atms) == 3:
            dih_a, dih_b, dih_c = adj_atms
            dih_ang = automol.geom.dihedral_angle(
                geo, dih_a, rad_atm, dih_b, dih_c)
            planarity = min([
                abs(dih_ang + x - numpy.pi)
                for x in [-2*numpy.pi, 0, 2*numpy.pi, 4*numpy.pi]])
            offset = [
                x for x in [-2*numpy.pi, 0, 2*numpy.pi, 4*numpy.pi]
                if abs(abs(dih_ang + x - numpy.pi) - planarity) < .01]
            start_inversion = dih_ang + offset - numpy.pi
            print(
                'checking for umbrella inversion around atom',
                rad_atm, start_inversion)
        #dont_dbl = ()
        #for rotor in rotors:
        #    for torsion in rotor:
        #        if umb_inv:
        #            break
        #        if any([ax_atm in dont_dbl for ax_atm in torsion.axis]):
        #            continue
        #        if any([rad_atm in torsion.axis for rad_atm in rad_atms]):
        #            for scan_geo in torsion.scan_geos.values():
        #                if scan_geo is None:
        #                    print('missing geo')
        #                    continue
        #                dih_ang = automol.geom.dihedral_angle(
        #                    scan_geo, dih_a, rad_atm, dih_b, dih_c)
        #                planarity = min([
        #                    abs(dih_ang + x - numpy.pi)
        #                    for x in [-2*numpy.pi, 0, 2*numpy.pi, 4*numpy.pi]])
        #                offset = [
        #                    x for x in [-2*numpy.pi, 0, 2*numpy.pi, 4*numpy.pi]
        #                    if abs(abs(dih_ang + x - numpy.pi) - planarity) < .01]
        #                inversion = dih_ang + offset - numpy.pi
        #                print(
        #                    'checking umbrella inversion d',
        #                    start_inversion, inversion)
        #                if start_inversion * inversion < 0:
        #                    #print(
        #                    #    'umbrella inversion occured',
        #                    #    start_inversion, inversion)
        #                    umb_inv = True
        #                    # break
        #            #if automol.pot.is_symmetric(torsion.pot):
        #            #    dont_dbl += tuple([rad_atm for rad_atm in torsion.axis if rad_atm in rad_atms])
        #            #    umb_fact *= 2
        #            #    ioprinter.info_message(
        #            #        'Umbrella mode identified for torsion about', torsion.axis)
    return 2 if planarity < .29 and planarity > .09 else 1
