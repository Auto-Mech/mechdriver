"""
  Therm Calculations
"""

import sys
import os
import math
import multiprocessing
import random
import automol.inchi
import automol.geom
from phydat import phycon
from mechanalyzer.inf import rxn as rinfo
from mechlib import filesys
from mechlib.filesys import reaction_fs
from mechlib.amech_io import printer as ioprinter
from mechroutines.pf.models.ene import read_energy
from mechroutines.pf.thermo import heatform


# FUNCTIONS TO PREPARE THE LIST OF REFERENCE SPECIES NEEDED FOR THERM CALCS #
REF_CALLS = {
    "basic": "get_basic",
    "cbh0": "get_cbhzed",
    "cbh1": "get_cbhone",
    "cbh1_0": "get_cbhone",
    "cbh1_1": "get_cbhone",
    "cbh2": "get_cbhtwo",
    "cbh2_0": "get_cbhtwo",
    "cbh2_1": "get_cbhtwo",
    "cbh2_2": "get_cbhtwo",
    "cbh3": "get_cbhthree"
}

TS_REF_CALLS = {
    "basic": "get_basic_ts",
    "cbh0": "get_cbhzed_ts",
    "cbh1": "get_cbhone_ts",
    "cbh1_0": "get_cbhzed_ts",
    "cbh1_1": "get_cbhone_ts",
    "cbh2": "get_cbhzed_ts",
    "cbh2_0": "get_cbhzed_ts",
    "cbh2_1": "get_cbhone_ts",
    "cbh3": "get_cbhone_ts"
}

IMPLEMENTED_CBH_TS_CLASSES = [
    'hydrogen abstraction high',
    # 'hydrogen migration',
    'beta scission',
    'elimination high',
    'radical radical hydrogen abstraction high',
    'addition high'
]


def prepare_refs(
        ref_scheme, spc_dct, spc_queue, run_prefix, save_prefix,
        repeats=False, parallel=False, zrxn=None):
    """ add refs to species list as necessary
    """
    spc_names = [spc[0] for spc in spc_queue]

    if parallel:
        nproc_avail = len(os.sched_getaffinity(0)) - 1

        num_spc = len(spc_names)
        spc_per_proc = math.floor(num_spc / nproc_avail)

        queue = multiprocessing.Queue()
        procs = []
        random.shuffle(spc_names)
        for proc_n in range(nproc_avail):
            spc_start = proc_n*spc_per_proc
            if proc_n == nproc_avail - 1:
                spc_end = num_spc
            else:
                spc_end = (proc_n+1)*spc_per_proc

            spc_lst = spc_names[spc_start:spc_end]

            proc = multiprocessing.Process(
                target=_prepare_refs,
                args=(queue, ref_scheme, spc_dct, spc_lst,
                      run_prefix, save_prefix,
                      repeats, parallel, zrxn))
            procs.append(proc)
            proc.start()

        basis_dct = {}
        unique_refs_dct = {}
        for _ in procs:
            bas_dct, unq_dct = queue.get()
            basis_dct.update(bas_dct)
            bas_ichs = [
                unique_refs_dct[spc]['inchi']
                if 'inchi' in unique_refs_dct[spc]
                else unique_refs_dct['reacs']
                for spc in unique_refs_dct]
            for spc in unq_dct:
                new_ich = (
                    unq_dct[spc]['inchi']
                    if 'inchi' in unq_dct[spc] else unq_dct[spc]['reacs'])
                if new_ich not in bas_ichs:
                    cnt = len(list(unique_refs_dct.keys())) + 1
                    if isinstance(new_ich, str):
                        ref_name = 'REF_{}'.format(cnt)
                        unique_refs_dct[ref_name] = unq_dct[spc]
                    else:
                        ref_name = 'TS_REF_{}'.format(cnt)
                        unique_refs_dct[ref_name] = unq_dct[spc]
        for proc in procs:
            proc.join()
    else:
        basis_dct, unique_refs_dct = _prepare_refs(
            None, ref_scheme, spc_dct, spc_names,
            run_prefix, save_prefix,
            repeats=repeats, parallel=parallel,
            zrxn=zrxn)

    return basis_dct, unique_refs_dct


def _prepare_refs(queue, ref_scheme, spc_dct, spc_names,
                  run_prefix, save_prefix,
                  repeats=False, parallel=False, zrxn=None):
    """ Prepare references
    """

    ioprinter.info_message(
        'Processor {} will prepare species: {}'.format(
            os.getpid(), ', '.join(spc_names)))
    spc_ichs = [spc_dct[spc]['inchi'] for spc in spc_names]
    dct_ichs = [spc_dct[spc]['inchi'] for spc in spc_dct.keys()
                if spc != 'global' and 'ts' not in spc]

    # Determine the function to be used to get the thermochemistry ref species
    if ref_scheme in REF_CALLS:
        get_ref_fxn = getattr(heatform, REF_CALLS[ref_scheme])
    if ref_scheme in TS_REF_CALLS:
        get_ts_ref_fxn = getattr(heatform, TS_REF_CALLS[ref_scheme])

    # Print the message
    msg = '\nDetermining reference molecules for scheme: {}'.format(ref_scheme)
    msg += '\n'

    basis_dct = {}
    unique_refs_dct = {}
    # ioprinter.info_message('spc dct: ', spc_dct.keys())
    # Determine the reference species, list of inchis
    for spc_name, spc_ich in zip(spc_names, spc_ichs):
        msg += '\nDetermining basis for species: {}'.format(spc_name)
        if zrxn is not None:
            rxnclass = automol.reac.reaction_class(zrxn)
            if (rxnclass in IMPLEMENTED_CBH_TS_CLASSES and
               'basic' not in ref_scheme):
                spc_basis, coeff_basis = get_ts_ref_fxn(zrxn)
            else:
                # Use a basic scheme
                spc_basis = []
                coeff_basis = []
                ts_ref_scheme = ref_scheme
                if '_' in ts_ref_scheme:
                    ts_ref_scheme = 'cbh' + ref_scheme.split('_')[1]
                for spc_i in spc_dct[spc_name]['reacs']:
                    bas_dct_i, _ = prepare_refs(
                        ts_ref_scheme, spc_dct, [[spc_i, None]],
                        run_prefix, save_prefix)
                    spc_bas_i, coeff_bas_i = bas_dct_i[spc_i]
                    for bas_i, c_bas_i in zip(spc_bas_i, coeff_bas_i):
                        if bas_i not in spc_basis:
                            spc_basis.append(bas_i)
                            coeff_basis.append(c_bas_i)
                        else:
                            for j, bas_j in enumerate(spc_basis):
                                if bas_i == bas_j:
                                    coeff_basis[j] += c_bas_i
        else:
            spc_basis, coeff_basis = get_ref_fxn(spc_ich)
        for i,  _ in enumerate(spc_basis):
            if isinstance(spc_basis[i], str):
                spc_basis[i] = automol.inchi.add_stereo(spc_basis[i])

        msg += '\nInCHIs for basis set:'
        for base in spc_basis:
            msg += '\n  {}'.format(base)

        # Add to the dct containing info on the species basis
        basis_dct[spc_name] = (spc_basis, coeff_basis)

        # Add to the dct with reference dct if it is not in the spc dct
        for ref in spc_basis:
            print('ref test', ref)
            bas_ichs = [
                unique_refs_dct[spc]['inchi']
                if 'inchi' in unique_refs_dct[spc]
                else unique_refs_dct[spc]['reacs']
                for spc in unique_refs_dct]
            cnt = len(list(unique_refs_dct.keys())) + 1
            if isinstance(ref, str):
                if ((ref not in spc_ichs and ref not in dct_ichs)
                        or repeats) and ref not in bas_ichs:
                    ref_name = 'REF_{}'.format(cnt)
                    msg += (
                        '\nAdding reference species {}, InChI string:{}'
                    ).format(ref, ref_name)
                    unique_refs_dct[ref_name] = create_spec(ref)
            else:
                if _chk(ref, spc_ichs, dct_ichs, bas_ichs, repeats):
                    ref_name = 'TS_REF_{}'.format(cnt)
                    msg += (
                        '\nAdding reference species {}, InChI string:{}'
                    ).format(ref, ref_name)
                    unique_refs_dct[ref_name] = create_ts_spc(
                        ref, spc_dct, spc_dct[spc_name]['mult'],
                        run_prefix, save_prefix,
                        rxnclass)
    ioprinter.info_message(msg)

    ret = None
    if parallel:
        queue.put((basis_dct, unique_refs_dct))
    else:
        ret = (basis_dct, unique_refs_dct)

    return ret


def create_ts_spc(ref, spc_dct, mult, run_prefix, save_prefix, rxnclass):
    """ add a ts species to the species dictionary
    """

    # Obtain the Reaction InChIs, Charges, Mults
    reacs, prods = ref[0], ref[1]
    rxn_ichs = (
        tuple(automol.inchi.add_stereo(ich) for ich in reacs if ich),
        tuple(automol.inchi.add_stereo(ich) for ich in prods if ich)
    )

    rxn_muls, rxn_chgs = (), ()
    for rgts in (reacs, prods):
        rgt_muls, rgt_chgs = (), ()
        for rgt in rgts:
            found = False
            for name in spc_dct:
                if 'inchi' in spc_dct[name]:
                    if spc_dct[name]['inchi'] == rgt:
                        rgt_muls += (spc_dct[name]['mult'],)
                        rgt_chgs += (spc_dct[name]['charge'],)
                        found = True
                        break
            if not found:
                new_spc = create_spec(rgt)
                rgt_muls += (new_spc['mult'],)
                rgt_chgs += (new_spc['charge'],)
        rxn_muls += (rgt_muls,)
        rxn_chgs += (rgt_chgs,)

    rxn_info = rinfo.from_data(rxn_ichs, rxn_chgs, rxn_muls, mult)

    return {
        'reacs': list(reacs),
        'prods': list(prods),
        'charge': 0,
        'inchi': '',
        'class': rxnclass,
        'mult': mult,
        'rxn_info': rxn_info,
        'ts_locs': (0,),
        'rxn_fs': reaction_fs(run_prefix, save_prefix, rinfo.sort(rxn_info))
    }


def create_spec(ich, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
    """ add a species to the species dictionary
    """
    rad = automol.formula.electron_count(automol.inchi.formula(ich)) % 2
    mult = 1 if not rad else 2

    return {
        'inchi': ich,
        'inchikey': automol.inchi.inchi_key(ich),
        'sens': 0.0,
        'charge': charge,
        'mult': mult,
        'mc_nsamp': mc_nsamp,
        'hind_inc': hind_inc * phycon.DEG2RAD
    }


def is_scheme(scheme):
    """ Return Boolean val if there is a scheme
    """
    return bool(scheme in REF_CALLS)


# FUNCTIONS TO CALCULATE ENERGIES FOR THERMOCHEMICAL PARAMETERS #
def basis_energy(spc_name, spc_basis, uni_refs_dct, spc_dct,
                 spc_model_dct_i, run_prefix, save_prefix,
                 read_species=True):
    """ Return the electronic + zero point energies for a set of species.
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

    basis_dct, uniref_dct = prepare_refs(
        ref_scheme, spc_dct, [[spc_name, None]],
        run_prefix, save_prefix,
        zrxn=zrxn)

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
    hf0k = heatform.calc_hform_0k(
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
                basis_dct_trs, uniref_dct_trs = prepare_refs(
                    ts_ref_scheme, spc_dct, [[spc_name, None]],
                    run_prefix, save_prefix,
                    zrxn=zrxn)
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
                hf0k_trs = heatform.calc_hform_0k(
                    ene_spc_trs, ene_basis_trs, spc_basis_trs,
                    coeff_basis_trs, ref_set=ref_enes)
            else:
                hf0k_trs = hf0k
        else:
            hf0k_trs = 0.0
    else:
        hf0k_trs = None
    ioprinter.info_message('ABS Energy  (kJ/mol): ', ene_chnlvl)
    ioprinter.info_message('Hf0K Energy (kJ/mol): ', hf0k * phycon.KCAL2KJ)

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


def _chk(ref, spc_ichs, dct_ichs, bas_ichs, repeats):
    """ a """

    ini = (
        ((ref not in spc_ichs and ref not in dct_ichs) or repeats) and
        (ref not in bas_ichs)
    )
    rref = ref[::-1]
    sec = (
        ((rref not in spc_ichs and rref not in dct_ichs) or repeats) and
        (ref not in bas_ichs[::-1])
    )

    return ini or sec
