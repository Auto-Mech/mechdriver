""" specifier mappings for naming directories in a DataSeries
"""
import os
import string
import numbers
import elstruct
import automol
from autofile._safemode import safemode_is_on
from autofile.schema._util import (is_valid_inchi_multiplicity as
                                   _is_valid_inchi_multiplicity)
from autofile.schema._util import short_hash as _short_hash
from autofile.schema._util import (random_string_identifier as
                                   _random_string_identifier)
from autofile.schema._util import (is_random_string_identifier as
                                   _is_random_string_identifier)


# Specifier mappings for species-specific layers
def species_trunk():
    """ species trunk directory name
    """
    return 'SPC'


def species_leaf(ich, chg, mul):
    """ species leaf directory name
    """
    if safemode_is_on():
        assert automol.chi.is_standard_form(ich)
        assert automol.chi.is_complete(ich)

    assert isinstance(chg, numbers.Integral)
    assert isinstance(mul, numbers.Integral), (
        f'Multiplicity {mul} is not an integer'
    )

    assert _is_valid_inchi_multiplicity(ich, mul, chg), (
        f'inchi is {ich} multiplicity is {mul} and charge is {chg}'
    )

    ick = automol.chi.inchi_key(ich)
    chg_str = str(chg)
    mul_str = str(mul)

    dir_names = (automol.chi.formula_layer(ich),
                 automol.inchi_key.first_hash(ick),
                 chg_str,
                 mul_str,
                 automol.inchi_key.second_hash_with_extension(ick))
    return os.path.join(*dir_names)


# Specifier mappings for reaction-specific layers
def reaction_trunk():
    """ reaction trunk directory name
    """
    return 'RXN'


def reaction_leaf(rxn_ichs, rxn_chgs, rxn_muls, ts_mul):
    """ reaction leaf directory name
    """
    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    if safemode_is_on():
        assert ((rxn_ichs, rxn_chgs, rxn_muls) ==
                sort_together(rxn_ichs, rxn_chgs, rxn_muls))
    ichs1, ichs2 = rxn_ichs
    chgs1, chgs2 = rxn_chgs
    muls1, muls2 = rxn_muls
    return os.path.join(_reactant_leaf(ichs1, chgs1, muls1),
                        _reactant_leaf(ichs2, chgs2, muls2),
                        str(ts_mul))


def reaction_is_reversed(rxn_ichs, rxn_chgs, rxn_muls):
    """ sort inchis, chgs, and muliplicities together
    """

    assert len(rxn_ichs) == len(rxn_chgs) == len(rxn_muls) == 2

    ichs1, ichs2 = rxn_ichs
    chgs1, chgs2 = rxn_chgs
    muls1, muls2 = rxn_muls

    ichs1, chgs1, muls1 = _sort_together(ichs1, chgs1, muls1)
    ichs2, chgs2, muls2 = _sort_together(ichs2, chgs2, muls2)

    return (_sortable_representation(ichs1, chgs1, muls1) >
            _sortable_representation(ichs2, chgs2, muls2))


def sort_together(rxn_ichs, rxn_chgs, rxn_muls):
    """ sort inchis, chgs, and muliplicities together
    """

    assert len(rxn_ichs) == len(rxn_chgs) == len(rxn_muls) == 2

    ichs1, ichs2 = rxn_ichs
    chgs1, chgs2 = rxn_chgs
    muls1, muls2 = rxn_muls

    ichs1, chgs1, muls1 = _sort_together(ichs1, chgs1, muls1)
    ichs2, chgs2, muls2 = _sort_together(ichs2, chgs2, muls2)

    if reaction_is_reversed(rxn_ichs, rxn_chgs, rxn_muls):
        ichs1, ichs2 = ichs2, ichs1
        chgs1, chgs2 = chgs2, chgs1
        muls1, muls2 = muls2, muls1

    return ((ichs1, ichs2), (chgs1, chgs2), (muls1, muls2))


def _sort_together(ichs, chgs, muls):
    idxs = automol.chi.argsort(ichs)
    ichs = tuple(ichs[idx] for idx in idxs)
    chgs = tuple(chgs[idx] for idx in idxs)
    muls = tuple(muls[idx] for idx in idxs)
    return (ichs, chgs, muls)


def _sortable_representation(ichs, chgs, muls):
    idxs = automol.chi.argsort(ichs)
    ichs = tuple(ichs[idx] for idx in idxs)
    return (len(ichs), ichs, chgs, muls)


def _reactant_leaf(ichs, chgs, muls):
    """ reactant leaf directory name
    """
    if safemode_is_on():
        for ich in ichs:
            assert automol.chi.is_standard_form(ich), (
                f'{ich} not standard form')
            assert automol.chi.is_complete(ich), (
                f'{ich} not complete')
        assert tuple(ichs) == automol.chi.sorted_(ichs)

    assert len(ichs) == len(chgs) == len(muls)
    assert all(isinstance(chg, numbers.Integral) for chg in chgs)
    assert all(isinstance(mul, numbers.Integral) for mul in muls)
    assert all(_is_valid_inchi_multiplicity(ich, mul)
               for ich, mul in zip(ichs, muls))

    ich = automol.chi.standard_form(automol.chi.join(ichs))
    ick = automol.chi.inchi_key(ich)
    chg_str = '_'.join(map(str, chgs))
    mul_str = '_'.join(map(str, muls))

    dir_names = (automol.chi.formula_layer(ich),
                 automol.inchi_key.first_hash(ick),
                 chg_str,
                 mul_str,
                 automol.inchi_key.second_hash_with_extension(ick))
    return os.path.join(*dir_names)


def transition_state_trunk():
    """ transition state trunk directory name
    """
    return 'TS'


def transition_state_leaf(num):
    """ transition state leaf directory name
    """
    assert isinstance(num, numbers.Integral) and 0 <= num <= 99, (
        f'Num {num} must be integer between 0 and 99'
    )
    return f'{int(num):02d}'


# Specifier mappings for layers used by both species and reaction file systems
def theory_leaf(method, basis, orb_type):
    """ theory leaf directory name

    This need not be tied to elstruct -- just take out the name checks.

    :param method: the name of the electronic structure method
    :type method: str
    :param basis: the atomic orbital basis set
    :type basis: str
    :param orb_type: 'R' indicates restricted orbitals, 'U' indicates
        unrestricted orbitals
    :type orb_type: str
    """

    # Pull apart the method
    core_method, pfxs = elstruct.Method.evaluate_method_type(method)

    # Build a hash of the prefixes
    if pfxs:
        ord_pfx_str = ''.join(list(pfxs))
        hashed_pfx = _short_hash(ord_pfx_str.lower()) + '-'
    else:
        hashed_pfx = ''

    assert elstruct.Method.contains(core_method)
    assert elstruct.Basis.contains(basis)

    if orb_type in ('R', 'U'):
        orb_hash = orb_type
    else:
        assert (
            isinstance(orb_type, list) and len(orb_type) == 2 and
            all(isinstance(x, int) for x in orb_type)
        ), (
            f"orb_type can only be 'R', 'U', or [int, int]\n"
            f"\tvalue received: {orb_type}")

        orb_hash = '@' + _short_hash(orb_type).upper()

    dir_name = ''.join([hashed_pfx,
                        _short_hash(core_method.lower()),
                        _short_hash(None if basis is None else basis.lower()),
                        orb_hash])
    return dir_name


def conformer_trunk():
    """ conformer trunk directory name
    """
    return 'CONFS'


def conformer_branch(rid):
    """ ring conformer leaf directory name
    """
    assert rid[0] == 'r', (
        f'rid {rid} does not start with r'
    )
    assert _is_random_string_identifier(rid[1:])
    return rid


def conformer_leaf(cid):
    """ torsion conformer leaf directory name
    """
    assert cid[0] == 'c', (
        f'cid {cid} does not start with c'
    )
    assert _is_random_string_identifier(cid[1:])

    return cid


def generate_new_conformer_id():
    """ generate a new conformer identifier
    """
    return 'c'+_random_string_identifier()


def generate_new_ring_id():
    """ generate a new conformer identifier
    """
    return 'r'+_random_string_identifier()


def single_point_trunk():
    """ single point trunk directory name
    """
    return 'SP'


def high_spin_trunk():
    """ high spin, single point trunk directory name
    """
    return 'HS'


def symmetry_trunk():
    """ symmetric-conformer trunk directory name
    """
    return 'SYM'


def zmatrix_trunk():
    """ zmatrix trunk directory name
    """
    return 'Z'


def zmatrix_leaf(num):
    """ zmatrix leaf directory name
    """
    assert isinstance(num, numbers.Integral) and 0 <= num <= 99, (
        f'Num {num} must be integer between 0 and 99'
    )
    return f'{int(num):02d}'


def scan_trunk():
    """ scan trunk directory name
    """
    return 'SCANS'


def scan_branch(coo_names):
    """ scan branch directory name
    """
    return '_'.join(sorted(coo_names))


def scan_leaf(coo_vals):
    """ scan leaf directory name
    """
    return '_'.join((f'{val:.2f}' for val in coo_vals))


def cscan_trunk():
    """ constrained scan trunk directory name
    """
    return 'CSCANS'


def cscan_branch1(cons_coo_val_dct):
    """ constrained scan branch 1 directory name

    :param cons_coo_val_dct: a dictionary of the constraint values, keyed by
        coordinate name
    :type cons_coo_val_dct: dict
    """
    coo_typ_sort_dct = {'R': 0, 'A': 1, 'D': 2}

    def _sort_key(coo_name):
        typ = coo_name[0]
        num = int(coo_name[1:])
        return coo_typ_sort_dct[typ], num

    # Sort the coordinate names by type and number
    cons_coo_names = sorted(cons_coo_val_dct.keys(), key=_sort_key)

    # Get the coordinate values in sorted order and round them
    cons_coo_vals = [cons_coo_val_dct[name] for name in cons_coo_names]
    cons_coo_vals = [float(round(val, 2)) for val in cons_coo_vals]

    cons_coo_val_dct = dict(zip(cons_coo_names, cons_coo_vals))
    return _short_hash(cons_coo_val_dct)


def cscan_branch2(coo_names):
    """ constrained scan branch 2 directory name
    """
    return '_'.join(sorted(coo_names))


def cscan_leaf(coo_vals):
    """ constrained scan leaf directory name
    """
    return '_'.join((f'{val:.2f}' for val in coo_vals))


def tau_trunk():
    """ tau trunk directory name
    """
    return 'TAU'


def tau_leaf(tid):
    """ tau leaf directory name
    """
    assert tid[0] == 't'
    assert _is_random_string_identifier(tid[1:])
    return tid


def generate_new_tau_id():
    """ generate a new conformer identifier
    """
    return 't'+_random_string_identifier()


def energy_transfer_trunk():
    """ energy_transfer trunk directory name
    """
    return 'ETRANS'


def vrctst_trunk():
    """ vrctst trunk directory name
    """
    return 'VRC'


def vrctst_leaf(num):
    """ vrctst leaf directory name
    """
    assert isinstance(num, numbers.Integral) and 0 <= num <= 99, (
        f'Num {num} must be integer between 0 and 99'
    )
    return f'{int(num):02d}'


# Specifier mappings specific to the run file system
def run_trunk():
    """ run trunk directory name
    """
    return 'RUN'


def run_leaf(job):
    """ run leaf directory name
    """
    assert elstruct.Job.contains(job)
    dir_name = job[:4].upper()
    return dir_name


def subrun_leaf(macro_idx, micro_idx):
    """ run leaf directory name
    """
    assert isinstance(macro_idx, numbers.Integral)
    assert isinstance(micro_idx, numbers.Integral)
    assert macro_idx < 26  # for now -- if needed we can add AA, AB, etc.
    macro_str = string.ascii_uppercase[macro_idx]
    micro_str = f'{micro_idx:0>2d}'
    return ''.join([macro_str, micro_str])


def build_trunk(head):
    """ build trunk directory name
    """
    assert isinstance(head, str)
    return head.upper()[:4]


def build_branch(_string):
    """ build branch directory name
    """
    assert isinstance(_string, str)
    return _string.upper()


def build_leaf(val):
    """ build leaf directory name
    """
    return str(val)


def get_next_build_number(num):
    """ determine the next build number
    """
    return int(num % 10)
