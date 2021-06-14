""" Handle symmetry factor stuff
"""

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

    sym_factor = spc_dct_i.get('sym_factor')
    if sym_factor is not None:
        ioprinter.info_message(
            ' - Reading symmetry number input by user:', sym_factor)
    else:

        zrxn = spc_dct_i.get('zrxn', None)
        if zrxn is not None:
            grxn = automol.reac.relabel_for_geometry(zrxn)
        else:
            grxn = None

        sym_model = spc_mod_dct_i['symm']['mod']

        # Obtain geometry, energy, and symmetry filesystem
        [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['symm']
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)

        # Obtain the external symssetry number
        ext_symm = automol.geom.external_symmetry_factor(geo)

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':

            # Set up the symmetry filesystem, read symmetrically similar geos
            # includes minimum geo
            sym_fs = fs.symmetry(cnf_path)
            symm_geos = [geo]
            symm_geos += [sym_fs[-1].file.geometry.read(locs)
                          for locs in sym_fs[-1].existing()]

            # Obtain the internal symmetry number and end group factors
            if rotors:
                ioprinter.info_message(
                    ' - Determining internal sym number ',
                    'using sampling routine.')
                int_symm, endgrp = automol.geom.internal_symm_from_sampling(
                    symm_geos, rotors, grxn=grxn, zma=zma)
            else:
                ioprinter.info_message(' - No torsions, internal sym is 1.0')
                int_symm, endgrp = 1.0, 1.0

        else:
            ioprinter.info_message(
                'No symmetry model requested, ',
                'setting internal sym factor to 1.0')
            int_symm, endgrp = 1.0, 1.0

        # Obtain overall number, reduced as needed
        int_symm = automol.geom.reduce_internal_symm(
            int_symm, ext_symm, endgrp, geo)

        sym_factor = ext_symm * int_symm

        sym_factor = automol.geom.rotor_reduced_symm_factor(
            sym_factor, automol.rotor.symmetries(rotors, flat=True))

    return sym_factor
