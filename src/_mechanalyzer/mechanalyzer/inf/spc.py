""" Constructs, Reads, and Builds Species Info Objects
"""

from mechanalyzer import par


SPC_PROPS = [
    par.SPC.INCHI,
    par.SPC.CHARGE,
    par.SPC.MULT
]

CANON_SPC_PROPS = [
    par.SPC.CANON_ENANT_ICH,
    par.SPC.CHARGE,
    par.SPC.MULT
]

def from_data(inchi, charge, mult):
    """ Construct a species info object from the constituent peieces of data.

        :param inchi: InChI string for species
        :type inchi: str
        :param charge: electric charge for species
        :type charge: int
        :param mult: spin multiplicity for species
        :type mult: int
        :rtype: tuple(str, int, int)
    """

    # assert to put in a data for a real species
    return (inchi, charge, mult)


def from_dct(dct, canonical=False):
    """ Construct a species info object by reading the constituent data
        from a dictionary containing the species data.

        :param dct: information dictionary for species
        :type dct: dict[str: obj]
        :rtype: tuple(str, int, int)
    """
    spc_props = SPC_PROPS
    if canonical:
        spc_props = CANON_SPC_PROPS
    assert set(spc_props) <= set(dct), (
        f'Properties {" ".join(SPC_PROPS)} not in dict'
    )

    inf_obj = tuple()
    for prop in spc_props:
        inf_obj += (dct[prop],)

    return inf_obj


def value(inf_obj, val):
    """ Obtain a desired value from a species info object.

        :param inf_obj: species info object
        :type inf_obj: mechanalyzer.inf.spc object
        :param val: value to obtain from info object
        :type val: str
        :rtype: str/int
    """

    assert val in SPC_PROPS, f'Desired value {val} not in spc info object'

    return inf_obj[SPC_PROPS.index(val)]


def combine(inf_obj1, inf_obj2):
    """ Created a spc info for two species where charge and mult
        have been combined
    """

    ich = (inf_obj1[0], inf_obj2[0])
    chg = inf_obj1[1] + inf_obj2[1]
    mult = max([inf_obj1[2], inf_obj2[2]])  # Taking max, assuming a complex

    return (ich, chg, mult)
