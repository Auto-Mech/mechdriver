""" elstruct.writer._mrcc2018 parameters
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
    lambda osp: 'scfmaxit={}'.format(*option.values(osp)),
    option.name(Option.Opt.MAXITER_):
    lambda osp: 'optmaxit'.format(*option.values(osp)),
}
