""" directory naming functions
"""
import os
import elstruct
import automol
from autofile.system._util import (is_valid_stereo_inchi as
                                   _is_valid_stereo_inchi)
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


def species_leaf(ich, mult):
    """ species leaf directory name
    """
    assert _is_valid_stereo_inchi(ich)
    assert _is_valid_inchi_multiplicity(ich, mult)
    ich_key = automol.inchi.inchi_key(ich)
    assert automol.inchi.key.is_standard_neutral(ich_key)

    mult_str = '{:d}'.format(mult)
    dir_names = (automol.inchi.formula_layer(ich),
                 automol.inchi.key.first_hash(ich_key),
                 mult_str,
                 automol.inchi.key.second_hash(ich_key),)
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
