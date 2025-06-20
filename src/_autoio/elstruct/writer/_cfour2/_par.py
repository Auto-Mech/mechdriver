""" elstruct.writer._cfour2 parameters
"""

from elstruct import Option
from elstruct import option
from elstruct.par import Reference


REF_DCT = {
    Reference.RHF: 'rhf',
    Reference.UHF: 'uhf'
}


OPTION_EVAL_DCT = {
    option.name(Option.Scf.MAXITER_):
    lambda osp: 'SCF_MAXCYC={}'.format(*option.values(osp)),
    option.name(Option.Opt.MAXITER_):
    lambda osp: 'GEO_MAXCYC={}'.format(*option.values(osp)),
}
