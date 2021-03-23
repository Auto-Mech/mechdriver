""" Handle symmetry factor stuff
"""

import automol
from autofile import fs
from mechlib.amech_io import printer as ioprinter


def symmetry_factor(pf_filesystems, pf_models, spc_dct_i, rotors, grxn=None):
    """ Calculate the symmetry factor for a species
        Note: ignoring for saddle pts the possibility that two configurations
        differ only in their torsional values.
        As a result, the symmetry factor is a lower bound of the true value
    """

    if 'sym_factor' in spc_dct_i:

        sym_factor = spc_dct_i['sym_factor']
        ioprinter.info_message(' - Reading symmetry number input by user:', sym_factor)

    else:

        zrxn = spc_dct_i.get('zrxn', None)
        if zrxn is not None:
            grxn = automol.reac.relabel_for_geometry(zrxn)
        else:
            grxn = None

        sym_model = pf_models['sym']

        # Obtain geometry, energy, and symmetry filesystem
        [cnf_fs, cnf_path, min_cnf_locs, _, _] = pf_filesystems['sym']
        geo = cnf_fs[-1].file.geometry.read(min_cnf_locs)

        # Obtain the external symssetry number
        ext_sym = automol.geom.external_symmetry_factor(geo)

        # Obtain the internal symmetry number using some routine
        if sym_model == 'sampling':

            # Set up the symmetry filesystem
            sym_fs = fs.symmetry(cnf_path)
            sym_geos = [geo]
            sym_geos += [sym_fs[-1].file.geometry.read(locs)
                         for locs in sym_fs[-1].existing()]

            # Obtain the internal
            if rotors:
                ioprinter.info_message(
                    ' - Determining internal sym number ',
                    'using sampling routine.')
                int_sym = automol.geom.internal_symmetry_number_from_sample(
                    sym_geos, grxn=grxn)
            else:
                ioprinter.info_message(' - No torsions, internal sym is 1.0')
                int_sym = 1.0

        else:
            ioprinter.info_message(
                'No symmetry model requested, ',
                'setting internal sym factor to 1.0')
            int_sym = 1.0

        # Obtain overall number
        sym_factor = ext_sym * int_sym

        # Reduce sym factor using rotor symmetries
        tors_symms = automol.rotor.symmetries(rotors, flat=True)
        for symm in tors_symms:
            sym_factor /= symm

    return sym_factor