"""
  Therm Calculations
"""

import automol.inchi
import automol.geom
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
        cnt = 1
        for ref in spc_basis:
            if ref not in spc_ichs:
                msg += 'Adding reference species ref_{} to dct\n'.format(ref)
                ref_name = 'REF_{}'.format(cnt)
                unique_refs_dct[ref_name] = create_spec(ref)

    return basis_dct, unique_refs_dct, msg


def create_spec(ich, charge=0,
                mc_nsamp=(True, 3, 1, 3, 100, 12),
                hind_inc=30.):
    """ add a species to the species dictionary
    """
    spec = {}
    rad = automol.formula.electron_count(automol.inchi.formula_dct(ich)) % 2
    mult = 1 if not rad else 2
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
def basis_energy(spc_bas, spc_dct):
    """ Return the electronic + zero point energies for a set of species.
    """
    h_basis = []
    for ich in spc_bas:
        for key in spc_dct:
            if ich == spc_dct[key]['ich']:
                ene = spc_dct[key]['ene'] + spc_dct[key]['zpe']/phycon.EH2KCAL
                h_basis.append(ene)
                break
    return h_basis


def get_hf0k(spc, spc_dct, spc_bas, coeff, ref_set='ANL0'):
    """ Determine the 0 K heat of formation from the
        species dictionary and a set of references species.
    """
    spc_ene = spc_dct[spc]['ene'] + spc_dct[spc]['zpe']/phycon.EH2KCAL
    h_basis = basis_energy(spc_bas, spc_dct)
    print('hf0k test:', spc, spc_dct[spc]['ene'], spc_dct[spc]['zpe'], spc_ene,
          spc_bas, h_basis, coeff)

    # Get the 0 K heat of formation
    # ref_set should be a parameter for this routine
    h0form = heatform.calc_hform_0k(
        spc_ene, h_basis, spc_bas, coeff, ref_set=ref_set)
    return h0form
