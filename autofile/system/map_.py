""" directory naming functions
"""
import os
import string
import numbers
import elstruct
import automol
from autofile.system._util import (is_valid_inchi_multiplicity as
                                   _is_valid_inchi_multiplicity)
from autofile.system._util import short_hash as _short_hash
from autofile.system._util import (random_string_identifier as
                                   _random_string_identifier)
from autofile.system._util import (is_random_string_identifier as
                                   _is_random_string_identifier)


# species
def species_trunk():
    """ species trunk directory name
    """
    return 'SPC'


def species_leaf(ich, charge, mult):
    """ species leaf directory name
    """
    return _reactant_leaf([ich], [charge], [mult])


# reactions
def reaction_trunk():
    """ reaction trunk directory name
    """
    return 'RXN'


def reaction_leaf(rxn_ichs, rxn_chgs, rxn_muls, ts_mult):
    """ reaction leaf directory name
    """
    rxn_ichs = tuple(map(tuple, rxn_ichs))
    rxn_chgs = tuple(map(tuple, rxn_chgs))
    rxn_muls = tuple(map(tuple, rxn_muls))
    assert ((rxn_ichs, rxn_chgs, rxn_muls) ==
            sort_together(rxn_ichs, rxn_chgs, rxn_muls))
    ichs1, ichs2 = rxn_ichs
    charges1, charges2 = rxn_chgs
    mults1, mults2 = rxn_muls
    return os.path.join(_reactant_leaf(ichs1, charges1, mults1),
                        _reactant_leaf(ichs2, charges2, mults2),
                        str(ts_mult))


def reaction_direction(rxn_ichs, rxn_chgs, rxn_muls):
    """ sort inchis, chgs, and muliplicities together
    """

    assert len(rxn_ichs) == len(rxn_chgs) == len(rxn_muls) == 2

    ichs1, ichs2 = rxn_ichs
    chgs1, chgs2 = rxn_chgs
    muls1, muls2 = rxn_muls

    ichs1, chgs1, muls1 = _sort_together(ichs1, chgs1, muls1)
    ichs2, chgs2, muls2 = _sort_together(ichs2, chgs2, muls2)

    return (_sortable_representation(ichs1, chgs1, muls1) <
            _sortable_representation(ichs2, chgs1, muls1))


def sort_together(rxn_ichs, rxn_chgs, rxn_muls):
    """ sort inchis, chgs, and muliplicities together
    """

    assert len(rxn_ichs) == len(rxn_chgs) == len(rxn_muls) == 2

    ichs1, ichs2 = rxn_ichs
    chgs1, chgs2 = rxn_chgs
    muls1, muls2 = rxn_muls

    ichs1, chgs1, muls1 = _sort_together(ichs1, chgs1, muls1)
    ichs2, chgs2, muls2 = _sort_together(ichs2, chgs2, muls2)

    if (_sortable_representation(ichs1, chgs1, muls1) >
            _sortable_representation(ichs2, chgs1, muls1)):
        ichs1, ichs2 = ichs2, ichs1
        chgs1, chgs2 = chgs2, chgs1
        muls1, muls2 = muls2, muls1

    return ((ichs1, ichs2), (chgs1, chgs2), (muls1, muls2))


def _sort_together(ichs, chgs, muls):
    idxs = automol.inchi.argsort(ichs)
    ichs = tuple(ichs[idx] for idx in idxs)
    chgs = tuple(chgs[idx] for idx in idxs)
    muls = tuple(muls[idx] for idx in idxs)
    return (ichs, chgs, muls)


def _sortable_representation(ichs, chgs, muls):
    return (len(ichs), sorted(automol.inchi.argsort(ichs)), chgs, muls)


def _reactant_leaf(ichs, charges, mults):
    """ reactant leaf directory name
    """
    assert all(map(automol.inchi.is_standard_form, ichs))
    assert all(map(automol.inchi.is_complete, ichs))
    assert tuple(ichs) == automol.inchi.sorted_(ichs)
    assert len(ichs) == len(charges) == len(mults)
    assert all(isinstance(charge, numbers.Integral) for charge in charges)
    assert all(isinstance(mult, numbers.Integral) for mult in mults)
    assert all(_is_valid_inchi_multiplicity(ich, mult)
               for ich, mult in zip(ichs, mults))

    ich = automol.inchi.standard_form(automol.inchi.join(ichs))
    ick = automol.inchi.inchi_key(ich)
    charge_str = '_'.join(map(str, charges))
    mult_str = '_'.join(map(str, mults))

    dir_names = (automol.inchi.formula_sublayer(ich),
                 automol.inchi_key.first_hash(ick),
                 charge_str,
                 mult_str,
                 automol.inchi_key.second_hash_with_extension(ick))
    return os.path.join(*dir_names)


# theory
def theory_leaf(method, basis, orb_restricted):
    """ theory leaf directory name

    This need not be tied to elstruct -- just take out the name checks.
    Note that we are (no longer) checking the orbital restriction.
    """
    assert elstruct.Method.contains(method)
    assert elstruct.Basis.contains(basis)
    assert isinstance(orb_restricted, bool)

    ref_char = 'R' if orb_restricted else 'U'
    dir_name = ''.join([_short_hash(method),
                        _short_hash(basis),
                        ref_char])
    return dir_name


# run
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
    micro_str = '{:0>2d}'.format(micro_idx)
    return ''.join([macro_str, micro_str])


# conformer
def conformer_trunk():
    """ conformer trunk directory name
    """
    return 'CONFS'


def conformer_leaf(cid):
    """ conformer leaf directory name
    """
    assert _is_random_string_identifier(cid)
    return cid


def generate_new_conformer_id():
    """ generate a new conformer identifier
    """
    return _random_string_identifier()


# single point
def single_point_trunk():
    """ single point trunk directory name
    """
    return 'SP'


# scan
def scan_trunk():
    """ scan trunk directory name
    """
    return 'SCANS'


def scan_branch(tors_names):
    """ scan branch directory name
    """
    return '_'.join(sorted(tors_names))


def scan_leaf(grid_idxs):
    """ scan leaf directory name
    """
    return '_'.join(map('{:0>2d}'.format, grid_idxs))


# tau
def tau_trunk():
    """ tau trunk directory name
    """
    return 'TAU'


def tau_leaf(cid):
    """ tau leaf directory name
    """
    assert _is_random_string_identifier(cid)
    return cid


# builds (MESS, NASA Poly, etc.)
def build_trunk(head):
    """ build trunk directory name
    """
    assert isinstance(head, str)
    return head.upper()[:4]


def build_leaf(num):
    """ build leaf directory name
    """
    assert isinstance(num, numbers.Integral) and 0 <= num <= 9
    return str(num)


def get_next_build_number(num):
    """ determine the next build number
    """
    return int(num % 10)
