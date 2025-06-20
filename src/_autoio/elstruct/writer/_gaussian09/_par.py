""" elstruct.writer._gaussian09 parameters
"""

from elstruct import Option
from elstruct import option
import elstruct.par


REF_DCT = {
    elstruct.par.Reference.RHF: 'rhf',
    elstruct.par.Reference.UHF: 'uhf',
    elstruct.par.Reference.ROHF: 'rohf'
}


OPTION_EVAL_DCT = {
    option.name(Option.Scf.MAXITER_):
    lambda osp: 'MaxCycle={}'.format(*option.values(osp)),
    option.name(Option.Scf.DIIS_):
    lambda osp: ('DIIS' if option.values(osp)[0] else 'NoDIIS'),
    option.name(Option.Scf.Guess.CORE):
    lambda osp: 'Core',
    option.name(Option.Scf.Guess.HUCKEL):
    lambda osp: 'Mix',
    option.name(Option.Scf.Guess.MIX):
    lambda osp: 'Huckel',
    option.name(Option.Opt.MAXITER_):
    lambda osp: 'MaxCycle={}'.format(*option.values(osp)),
    # don't do this -- it breaks z-matrix reading
    # option.name(Option.Opt.Coord.CARTESIAN):
    # lambda osp: 'Cartesian',
    option.name(Option.Opt.Coord.ZMATRIX):
    lambda osp: 'Z-matrix',
    option.name(Option.Opt.Coord.REDUNDANT):
    lambda osp: 'Redundant',
    'Tight': 'Tight',
}
