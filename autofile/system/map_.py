""" directory naming functions
"""
import os
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


def reaction_leaf(ichs_pair, charges_pair, mults_pair):
    """ reaction leaf directory name
    """
    assert ((ichs_pair, charges_pair, mults_pair) ==
            sort_together(ichs_pair, charges_pair, mults_pair))
    ichs1, ichs2 = ichs_pair
    charges1, charges2 = charges_pair
    mults1, mults2 = mults_pair
    return os.path.join(_reactant_leaf(ichs1, charges1, mults1),
                        _reactant_leaf(ichs2, charges2, mults2))


def sort_together(ichs_pair, charges_pair, mults_pair):
    """ sort inchis, charges, and multiplicities together
    """

    def _sort_together(ichs, charges, mults):
        idxs = automol.inchi.argsort(ichs)
        ichs = tuple(ichs[idx] for idx in idxs)
        charges = tuple(charges[idx] for idx in idxs)
        mults = tuple(mults[idx] for idx in idxs)
        return (ichs, charges, mults)

    def _sortable_representation(ichs):
        return (len(ichs), sorted(automol.inchi.argsort(ichs)))

    assert len(ichs_pair) == len(charges_pair) == len(mults_pair) == 2

    ichs1, ichs2 = ichs_pair
    charges1, charges2 = charges_pair
    mults1, mults2 = mults_pair

    ichs1, charges1, mults1 = _sort_together(ichs1, charges1, mults1)
    ichs2, charges2, mults2 = _sort_together(ichs2, charges2, mults2)

    if _sortable_representation(ichs1) > _sortable_representation(ichs2):
        ichs1, ichs2 = ichs2, ichs1
        charges1, charges2 = charges2, charges1
        mults1, mults2 = mults2, mults1

    return ((ichs1, ichs2), (charges1, charges2), (mults1, mults2))


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
