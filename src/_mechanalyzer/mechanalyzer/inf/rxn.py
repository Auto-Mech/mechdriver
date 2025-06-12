""" Constructs, Reads, and Builds Reaction Info Objects
"""

from autofile.schema import sort_together
from mechanalyzer import par
from mechanalyzer.inf import spc


RXN_PROPS = [
    par.SPC.INCHI,
    par.SPC.CHARGE,
    par.SPC.MULT,
    par.SPC.TSMULT
]


def from_dct(reacs, prods, spc_dct, rxn_mul='low'):
    """ Construct a full reaction info object using the names of the
        reactants and products, which are used to read a species dictionary for
        the required physical information.

        :param reacs: names of the reactants
        :type reacs: tuple(str)
        :param prods: names of the products
        :type prods: tuple(str)
        :param spc_dct:
        :type spc_dct: dict[]
        :param rxn_mul: multiplicity of reaction to store in object
        :type rxn_mul: str
    """

    # Build the tuples of the reacs+prods infos
    rxn_ichs, rxn_chgs, rxn_muls = tuple(), tuple(), tuple()
    for side in (reacs, prods):
        ichs, chgs, muls = tuple(), tuple(), tuple()
        for rct in side:
            spc_info = spc.from_dct(spc_dct[rct])
            ichs += (spc.value(spc_info, par.SPC.INCHI),)
            chgs += (spc.value(spc_info, par.SPC.CHARGE),)
            muls += (spc.value(spc_info, par.SPC.MULT),)
        rxn_ichs += (ichs,)
        rxn_chgs += (chgs,)
        rxn_muls += (muls,)

    # Determine the multiplicity of the full reaction
    fake_rxn_info = (rxn_ichs, rxn_chgs, rxn_muls, ())
    ts_mul = ts_mult(fake_rxn_info, rxn_mul=rxn_mul)

    rxn_info = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return rxn_info


def from_data(rxn_ichs, rxn_chgs, rxn_muls, ts_mul):
    """ Constructs a full reaction info object using the constituent data.

        :param rxn_ichs: InChI strings for reactants and products
        :type rxn_ichs: tuple(tuple(str), tuple(str))
        :param rxn_chgs: electric charges of reactants and products
        :type rxn_chgs: tuple(tuple(int), tuple(int))
        :param rxn_muls: spin multiplicities of reactants and products
        :type rxn_muls: tuple(tuple(int), tuple(int))
        :param ts_mul: spin multiplicity of reaction transition state
        :type ts_mul: int
    """
    return (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)


def value(inf_obj, val):
    """ Obtain a desired value from a reaction info object.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :param val: value to obtain from info object
        :type val: str
        :rtype: str/int
    """

    assert val in RXN_PROPS, (
        f'Desired value {val} not in rxn info object'
    )

    return inf_obj[RXN_PROPS.index(val)]


def ts_info(inf_obj):
    """ Construct a species info object for the reaction transisiton state.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: mechanalyzer.inf.spc object
    """

    _chg = ts_chg(inf_obj)
    _mul = value(inf_obj, par.SPC.TSMULT)

    return ('', _chg, _mul)


def rgts_info(inf_obj):
    """ Construct species info objects for all of the species
        that make up the reactants and products of the reaction.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: tuple(tuple(mechanalyzer.inf.spc object))
    """

    _rgts_info = ()
    for rgt in ('reacs', 'prods'):
        _rgts_info += (rgt_info(inf_obj, rgt),)

    return _rgts_info


def rgt_info(inf_obj, rgt):
    """ Construct species info objects for all of the species
        that make up one side of the reaction

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: tuple(mechanalyzer.inf.spc object)
    """

    assert rgt in ('reacs', 'prods')

    rxn_ichs, rxn_chgs, rxn_muls, _ = inf_obj
    if rgt == 'reacs':
        rgt_ichs, rgt_chgs, rgt_muls = rxn_ichs[0], rxn_chgs[0], rxn_muls[0]
    else:
        rgt_ichs, rgt_chgs, rgt_muls = rxn_ichs[1], rxn_chgs[1], rxn_muls[1]

    _rgt_info = ()
    for rgt_ich, rgt_chg, rgt_mul in zip(rgt_ichs, rgt_chgs, rgt_muls):
        _rgt_info += (spc.from_data(rgt_ich, rgt_chg, rgt_mul),)

    return _rgt_info


def sort(inf_obj, scheme='autofile'):
    """ Construct a new reaction object where the reaction direction and
        constituent reactants and products have been sorted according
        to some requested scheme.

        The default sorting corresponds to the autofile scheme used
        for building filesystem paths.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :param scheme: sorting scheme
        :type scheme: str
        :rtype: mechanalyzer.inf.rxn object
    """

    rxn_ichs, rxn_chgs, rxn_muls, ts_mul = inf_obj

    if scheme == 'autofile':
        rxn_ichs, rxn_chgs, rxn_muls = sort_together(
            rxn_ichs, rxn_chgs, rxn_muls)

    sort_inf_obj = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return sort_inf_obj


def reverse(inf_obj):
    """ Construct a new reaction info object for the reverse reaction
        where all of the reactant and product info has been flipped.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: mechanalyzer.inf.rxn object
    """

    rxn_ichs, rxn_chgs, rxn_muls, ts_mul = inf_obj

    rxn_ichs = (rxn_ichs[1], rxn_ichs[0])
    rxn_chgs = (rxn_chgs[1], rxn_chgs[0])
    rxn_muls = (rxn_muls[1], rxn_muls[0])

    rev_inf_obj = (rxn_ichs, rxn_chgs, rxn_muls, ts_mul)

    return rev_inf_obj


def ts_chg(inf_obj):
    """ Evaulate the electric charge of the transition state, i.e.,
        net charge of the surface the reaction occurs on.

        Determines the value using the charges of the reactants and products.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: int
    """

    rxn_chgs = value(inf_obj, par.SPC.CHARGE)

    _chg = 0
    for rct_chg in rxn_chgs[0]:
        _chg += rct_chg

    return _chg


def ts_mult(inf_obj, rxn_mul='low'):
    """ Evaulate the multilicity of the transition state, i.e.,
        spin-state of the surface the reaction occurs on.

        Determines both the `low-spin` or `high-spin` value using
        the multiplicities of the reactants and products. Then
        returns the version that is requested.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: int
    """

    rxn_muls = value(inf_obj, par.SPC.MULT)

    nrcts, nprds = len(rxn_muls[0]), len(rxn_muls[1])

    # Set the multiplicities
    rct_spin_sum, prd_spin_sum = 0, 0
    rct_muls, prd_muls = [], []
    if nrcts == 1 and nprds == 1:
        _mul = max(rxn_muls[0][0], rxn_muls[1][0])
    else:
        for rct_mul in rxn_muls[0]:
            rct_spin_sum += (rct_mul - 1.)/2.
            rct_muls.append(rct_mul)
        for prd_mul in rxn_muls[1]:
            prd_spin_sum += (prd_mul - 1.)/2.
            prd_muls.append(prd_mul)

        if rxn_mul == 'low':
            _mul = min(rct_spin_sum, prd_spin_sum)
        else:
            _mul = max(rct_spin_sum, prd_spin_sum)
        _mul = int(round(2*_mul + 1))

    return _mul


def radrad(inf_obj):
    """ Determine if a reaction can be classified as a recombination of
        two radicals, either in the forward or reverse direction.

        :param inf_obj: reaction info object
        :type inf_obj: mechanalyzer.inf.rxn object
        :rtype: bool
    """

    muls = value(inf_obj, par.SPC.MULT)
    rct_muls = muls[0]
    if len(rct_muls) > 1:
        mul1, mul2 = rct_muls
        rad_rad = bool(mul1 > 1 and mul2 > 1)
    else:
        prd_muls = muls[1]
        if len(prd_muls) > 1:
            mul1, mul2 = prd_muls
            rad_rad = bool(mul1 > 1 and mul2 > 1)
        else:
            rad_rad = False

    return rad_rad
