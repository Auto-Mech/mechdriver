"""
  Therm Calculations
"""

import sys
import automol.inchi
import automol.geom
from routines.pf.messf.ene import get_fs_ene_zpe
from routines.pf.thermo import heatform
from lib.phydat import phycon


# FUNCTIONS TO PREPARE THE LIST OF REFERENCE SPECIES NEEDED FOR THERM CALCS #
REF_CALLS = {"basic": "get_basic",
             "cbh0": "get_cbhzed",
             "cbh1": "get_cbhone",
             "cbh2": "get_cbhtwo"}


def prepare_refs(ref_scheme, spc_dct, spc_queue):
    """ add refs to species list as necessary
    """
    # Get a lst of ichs corresponding to the spc queue
    spc_ichs = [spc_dct[spc[0]]['ich'] for spc in spc_queue]
    dct_ichs = [spc_dct[spc]['ich'] for spc in spc_dct.keys()
                if spc != 'global']

    # Determine the function to be used to get the thermochemistry ref species
    if ref_scheme in REF_CALLS:
        get_ref_fxn = getattr(heatform, REF_CALLS[ref_scheme])
    else:
        raise NotImplementedError

    # Print the message
    msg = 'Determining {} reference molecules for: \n'.format(ref_scheme)

    # Determine the reference species, list of inchis
    basis_dct = {}
    unique_refs_dct = {}
    for spc_name, spc_ich in zip(spc_queue, spc_ichs):
        # Determine basis set for each spc using the specified therm scheme
        spc_basis, coeff_basis = get_ref_fxn(spc_ich)
        msg += 'Species {} with basis {}\n'.format(
            spc_ich, ', '.join(spc_basis))

        # Add to the dct containing info on the species basis
        basis_dct[spc_name] = (spc_basis, coeff_basis)

        # Add to the dct with reference dct if it is not in the spc dct
        print('basis', spc_basis)
        print('ichs', spc_ichs)
        print('dct ichs', dct_ichs)
        cnt = 1
        for ref in spc_basis:
            if ref not in spc_ichs and ref not in dct_ichs:
                ref_name = 'REF_{}'.format(cnt)
                msg += 'Adding reference species {}, InChI string:{}\n'.format(
                    ref, ref_name)
                unique_refs_dct[ref_name] = create_spec(ref)
                cnt += 1

    return basis_dct, unique_refs_dct, msg


def create_spec(ich, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
    """ add a species to the species dictionary
    """
    spec = {}
    rad = automol.formula.electron_count(automol.inchi.formula_dct(ich)) % 2
    mult = 1 if not rad else 2
    print('ich', ich)
    print(automol.inchi.geometry(ich))
    spec['zmatrix'] = automol.geom.zmatrix(automol.inchi.geometry(ich))
    spec['ich'] = ich
    spec['chg'] = charge
    spec['mul'] = mult
    spec['mc_nsamp'] = mc_nsamp
    spec['hind_inc'] = hind_inc * phycon.DEG2RAD
    return spec


def is_scheme(scheme):
    """ Return Boolean val if there is a scheme
    """
    return bool(scheme in REF_CALLS)


# FUNCTIONS TO CALCULATE ENERGIES FOR THERMOCHEMICAL PARAMETERS #
def basis_energy(spc_bas, uni_refs_dct, spc_dct,
                 thy_dct, model_dct, model, save_prefix):
    """ Return the electronic + zero point energies for a set of species.
    """

    # Initialize ich name dct to noe
    ich_name_dct = {}
    for ich in spc_bas:
        ich_name_dct[ich] = None

    # Get names from the respective spc dcts
    for ich in spc_bas:
        for name in spc_dct:
            if name != 'global':
                if ich == spc_dct[name]['ich']:
                    ich_name_dct[ich] = name
        for name in uni_refs_dct:
            if ich == uni_refs_dct[name]['ich']:
                ich_name_dct[ich] = name

    # Check the ich_name_dct
    dct_incomplete = False
    for ich, name in ich_name_dct.items():
        if name is None:
            print('{} not given in species.csv file'.format(ich))
            dct_incomplete = True
    if dct_incomplete:
        print('*ERROR Job ending since basis species not specified')
        sys.exit()

    # Combine the two spc dcts together
    full_spc_dct = {**spc_dct, **uni_refs_dct}

    # Get the energies
    h_basis = []
    for ich, name in ich_name_dct.items():
        h_basis.append(
            get_fs_ene_zpe(
                full_spc_dct, name,
                thy_dct, model_dct, model,
                save_prefix, saddle=False,
                read_ene=True, read_zpe=True))

    # Check if all the energies found
    no_ene_cnt = 0
    for basis_ene, basis_name in zip(h_basis, ich_name_dct.values()):
        if basis_ene is None:
            print('No energy found for {}'.format(basis_name))
            no_ene_cnt += 1
    if no_ene_cnt > 1:
        print('*ERROR: Not all energies found for the basis species')
        sys.exit()

    return h_basis
