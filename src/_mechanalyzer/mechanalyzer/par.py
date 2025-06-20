"""
  Set various type parameters used in automol
"""


class SPC():
    """ Params for species
    """
    NAME = 'name'
    INCHI = 'inchi'
    CANON_ENANT_ICH = 'canon_enant_ich'
    INCHI_KEY = 'inchikey'
    SMILES = 'smiles'
    CHARGE = 'charge'
    MULT = 'mult'
    TSMULT = 'tsmult'
    SENS = 'sensitivity'
    ELECLVLS = 'electronic_levels'


class RXN():
    """ Params for reactions
    """
    REACS = 'reacs'
    PRODS = 'prods'


class THY():
    """ Params for thy
    """
    PROGRAM = 'program'
    METHOD = 'method'
    BASIS = 'basis'
    ORB_RESTRICT = 'orb_res'
