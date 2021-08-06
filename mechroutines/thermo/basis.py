"""
  Therm Calculations
"""

import sys
import automol.inchi
import automol.geom
from phydat import phycon
import thermfit
from mechlib import filesys
from mechlib.amech_io import printer as ioprinter
from mechroutines.models.ene import read_energy


# FUNCTIONS TO CALCULATE ENERGIES FOR THERMOCHEMICAL PARAMETERS #
def basis_energy(spc_name, spc_basis, uni_refs_dct, spc_dct,
                 spc_model_dct_i, run_prefix, save_prefix,
                 read_species=True):
    """ Reads the electronic and zero-point energies for a species and
        transition state and their constituent basis set.
    """

    # Initialize ich name dct to noe
    ich_name_dct = {}
    for ich in spc_basis:
        if isinstance(ich, str):
            ich_name_dct[ich] = None
        else:
            ich_name_dct[_ich_key_name(ich)] = None

    # Get names of the basis species from the respective spc dcts
    for ich in spc_basis:
        for name in spc_dct:
            if name != 'global' and 'ts' not in name:
                if ich == spc_dct[name]['inchi']:
                    ich_name_dct[ich] = name
            elif name != 'global':
                if 'reacs' in spc_dct[name]:
                    if _ich_in_rxn(ich, spc_dct[name]):
                        ich_name_dct[_ich_key_name(ich)] = name
        for name in uni_refs_dct:
            if 'TS' in name:
                ioprinter.info_message(uni_refs_dct[name]['reacs'])
                if _ich_in_rxn(ich, spc_dct[name]):
                    ich_name_dct[_ich_key_name(ich)] = name
            elif ich == uni_refs_dct[name]['inchi']:
                ich_name_dct[ich] = name

    # Check the ich_name_dct
    dct_incomplete = False
    for ich, name in ich_name_dct.items():
        if name is None:
            ioprinter.warning_message(
                '{} not given in species.csv file'.format(ich))
            dct_incomplete = True
    if dct_incomplete:
        ioprinter.error_message('Job ending since basis species not specified')
        sys.exit()

    # Get the species energy
    if read_species:
        ioprinter.info_message(
            'Calculating energy for species {}'.format(spc_name), newline=1)
        pf_filesystems = filesys.models.pf_filesys(
            spc_dct[spc_name], spc_model_dct_i,
            run_prefix, save_prefix, saddle='ts' in spc_name)
        h_spc = read_energy(
            spc_dct[spc_name], pf_filesystems,
            spc_model_dct_i, run_prefix,
            read_ene=True, read_zpe=True, saddle='ts' in spc_name)
        if h_spc is None:
            ioprinter.error_message('No energy found for {}'.format(spc_name))
            sys.exit()
    else:
        h_spc = None

    # Get the energies of the bases
    h_basis = []
    for ich, name in ich_name_dct.items():
        if name in spc_dct:
            spc_dct_i = spc_dct[name]
            prname = name
        elif name in uni_refs_dct:
            spc_dct_i = uni_refs_dct[name]
            prname = name
        if 'ts' in name or 'TS' in name:
            reacs, prods = ich.split('PRODS')
            reacs = reacs.replace('REACS', '')
            reacs = reacs.split('REAC')
            prods = prods.split('PROD')
            reac_lbl = 'r0'
            if len(reacs) > 1:
                reac_lbl += '+r1'
            prod_lbl = 'p0'
            if len(prods) > 1:
                prod_lbl += '+p1'
            ioprinter.info_message(
                'Basis Reaction: {}={} 1 1 1 '.format(reac_lbl, prod_lbl))
            for i, reac in enumerate(reacs):
                ioprinter.info_message(
                    'r{},{},{},{}'.format(
                        str(i), reac, automol.inchi.smiles(reac), '1'))
            for i, prod in enumerate(prods):
                ioprinter.info_message(
                    'p{},{},{},{}'.format(
                        str(i), prod, automol.inchi.smiles(prod), '1'))
        ioprinter.debug_message('bases energies test:', ich, name)
        pf_filesystems = filesys.models.pf_filesys(
            spc_dct_i, spc_model_dct_i,
            run_prefix, save_prefix, 'ts' in name or 'TS' in name)
        ioprinter.info_message(
            'Calculating energy for basis {}...'.format(prname), newline=1)
        h_basis.append(
            read_energy(
                spc_dct_i, pf_filesystems,
                spc_model_dct_i, run_prefix,
                read_ene=True, read_zpe=True,
                saddle='ts' in name or 'TS' in name
            )
        )

    # Check if all the energies found
    no_ene_cnt = 0
    for basis_ene, basis_name in zip(h_basis, ich_name_dct.values()):
        if basis_ene is None:
            ioprinter.warning_message(
                'No energy found for {}'.format(basis_name))
            no_ene_cnt += 1
    if no_ene_cnt > 1:
        ioprinter.error_message(
            'Not all energies found for the basis species')
        sys.exit()

    return h_spc, h_basis


def enthalpy_calculation(
        spc_dct, spc_name, ene_chnlvl,
        chn_basis_ene_dct, pes_mod_dct_i, spc_mod_dct_i,
        run_prefix, save_prefix,
        pforktp='ktp', zrxn=None):
    """ Calculate the Enthalpy.
    """

    ref_scheme = pes_mod_dct_i['therm_fit']['ref_scheme']
    ref_enes = pes_mod_dct_i['therm_fit']['ref_enes']

    basis_dct = thermfit.prepare_basis(
        ref_scheme, spc_dct, (spc_name,), zrxn=zrxn)
    uniref_dct = thermfit.unique_basis_species(basis_dct, spc_dct)

    # Get the basis info for the spc of interest
    spc_basis, coeff_basis = basis_dct[spc_name]

    ene_spc = ene_chnlvl
    ene_basis = []

    energy_missing = False
    for spc_basis_i in spc_basis:
        if not isinstance(spc_basis_i, str):
            basreacs, basprods = spc_basis_i
            spc_basis_i = ''.join(basreacs)
            spc_basis_i += ''.join(basprods)
        if spc_basis_i in chn_basis_ene_dct:
            ioprinter.debug_message(
                'Energy already found for basis species: ', spc_basis_i)
            ene_basis.append(chn_basis_ene_dct[spc_basis_i])
        else:
            ioprinter.debug_message(
                'Energy will be determined for basis species: ', spc_basis_i)
            energy_missing = True

    # Get the energies for the spc and its basis
    if energy_missing:
        _, ene_basis = basis_energy(
            spc_name, spc_basis, uniref_dct, spc_dct,
            spc_mod_dct_i,
            run_prefix, save_prefix, read_species=False)
        for spc_basis_i, ene_basis_i in zip(spc_basis, ene_basis):
            if not isinstance(spc_basis_i, str):
                basreacs, basprods = spc_basis_i
                spc_basis_i = ''.join(basreacs)
                spc_basis_i += ''.join(basprods)
            chn_basis_ene_dct[spc_basis_i] = ene_basis_i

    # Calculate and store the 0 K Enthalpy
    hf0k = thermfit.heatform.calc_hform_0k(
        ene_spc, ene_basis, spc_basis, coeff_basis, ref_set=ref_enes)

    if pforktp == 'ktp':
        if 'basic' in ref_scheme:
            ts_ref_scheme = 'basic'
        else:
            ts_ref_scheme = 'cbh0'
            if '_' in ref_scheme:
                ts_ref_scheme = 'cbh' + ref_scheme.split('_')[1]
        if zrxn is None:
            if ref_scheme != ts_ref_scheme:
                basis_dct_trs, uniref_dct_trs = thermfit.prepare_basis(
                    ts_ref_scheme, spc_dct, (spc_name,), zrxn=zrxn)
                uniref_dct_trs = thermfit.unique_basis_species(
                    basis_dct_trs, spc_dct)

                spc_basis_trs, coeff_basis_trs = basis_dct_trs[spc_name]
                ene_basis_trs = []
                energy_missing = False
                for spc_basis_i in spc_basis_trs:
                    if spc_basis_i in chn_basis_ene_dct:
                        ioprinter.info_message(
                            'Energy already found for basis species: ',
                            spc_basis_i)
                        ene_basis_trs.append(chn_basis_ene_dct[spc_basis_i])
                    else:
                        ioprinter.info_message(
                            'Energy will be determined for basis species: ',
                            spc_basis_i)
                        energy_missing = True
                if energy_missing:
                    _, ene_basis_trs = basis_energy(
                        spc_name, spc_basis_trs, uniref_dct_trs, spc_dct,
                        spc_mod_dct_i,
                        run_prefix, save_prefix)
                    for spc_basis_i, ene_basis_i in zip(
                            spc_basis_trs, ene_basis_trs):
                        chn_basis_ene_dct[spc_basis_i] = ene_basis_i
                ene_spc_trs = ene_chnlvl
                hf0k_trs = thermfit.heatform.calc_hform_0k(
                    ene_spc_trs, ene_basis_trs, spc_basis_trs,
                    coeff_basis_trs, ref_set=ref_enes)
            else:
                hf0k_trs = hf0k
        else:
            hf0k_trs = 0.0
    else:
        hf0k_trs = None
    ioprinter.info_message('ABS Energy  (hart): ', ene_chnlvl)
    ioprinter.info_message('Hf0K Energy (hart): ', hf0k * phycon.KCAL2KJ)

    return hf0k, hf0k_trs, chn_basis_ene_dct, basis_dct


# Helpers
def _ich_key_name(ich):
    """ Build a useful dictionary key of a joined ich
    """
    return ('REACS' + 'REAC'.join(ich[0]) +
            r'PRODS' + 'PROD'.join(ich[1]))


def _ich_in_rxn(ich, spc_dct_i):
    """ Check if in a reaction
    """

    reacs = spc_dct_i['reacs']
    reacs_rev = reacs[::-1]
    prods = spc_dct_i['prods']
    prods_rev = prods[::-1]

    return bool(
        list(ich[0]) in (reacs, reacs_rev) and
        list(ich[1]) in (prods, prods_rev)
    )
