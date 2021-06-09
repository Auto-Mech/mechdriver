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
        ioprinter.info_message(' - Reading symmetry number input by user:', sym_factor)
    else:

        zrxn = spc_dct_i.get('zrxn', None)
        if zrxn is not None:
            grxn = automol.reac.relabel_for_geometry(zrxn)
        else:
            grxn = None

        # if automol.geom.is_atom(geo):
        sym_model = spc_mod_dct_i['symm']['mod']

        # Obtain geometry, energy, and symmetry filesystem
        [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['symm']
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)

        # Obtain the external symssetry number
        ext_sym = automol.geom.external_symmetry_factor(geo)

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':

            # Set up the symmetry filesystem
            sym_fs = fs.symmetry(cnf_path)
            symm_geos = [geo]
            symm_geos += [sym_fs[-1].file.geometry.read(locs)
                         for locs in sym_fs[-1].existing()]

            # Obtain the internal
            if rotors:
                ioprinter.info_message(
                    ' - Determining internal sym number ',
                    'using sampling routine.')
                int_sym, end_group_factor = int_sym_num_from_sampling(
                    symm_geos, rotors, grxn=grxn, zma=zma)
            else:
                ioprinter.info_message(' - No torsions, internal sym is 1.0')
                int_sym = 1.0
                end_group_factor = 1.0

        else:
            ioprinter.info_message(
                'No symmetry model requested, ',
                'setting internal sym factor to 1.0')
            int_sym = 1.0
            end_group_factor = 1.0

        # Obtain overall number
        if ext_sym % 3 == 0 and end_group_factor > 1:
            if not automol.graph.is_branched(automol.geom.graph(geo)):
                int_sym = int_sym/3

        sym_factor = ext_sym * int_sym
        # Reduce sym factor using rotor symmetries
        sym_factor = tors_reduced_sym_factor(sym_factor, rotors)

    return sym_factor
