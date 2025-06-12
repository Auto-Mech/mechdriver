""" elstruct.writer._psi4 parameters
"""

import elstruct.option
import elstruct.par


REF_DCT = {
    elstruct.par.Reference.RHF: 'rhf',
    elstruct.par.Reference.UHF: 'uhf',
    elstruct.par.Reference.ROHF: 'rohf',
    elstruct.par.Reference.RKS: 'rks',
    elstruct.par.Reference.UKS: 'uks'
}

OPTION_EVAL_DCT = {
    elstruct.option.name(elstruct.par.Option.Mol.NOSYMM_):
    'symmetry c1',
    elstruct.option.name(elstruct.par.Option.Scf.MAXITER_):
    lambda osp: 'set scf maxiter {}'.format(*elstruct.option.values(osp)),
    elstruct.option.name(elstruct.par.Option.Scf.DIIS_):
    lambda osp: 'set scf diis {}'.format(*elstruct.option.values(osp)),
    elstruct.option.name(elstruct.par.Option.Scf.Guess.CORE):
    lambda osp: 'set scf guess core',
    elstruct.option.name(elstruct.par.Option.Scf.Guess.HUCKEL):
    lambda osp: 'set scf guess huckel',
    elstruct.option.name(elstruct.par.Option.Opt.MAXITER_):
    lambda osp: 'set geom_maxiter {}'.format(*elstruct.option.values(osp)),
    elstruct.option.name(elstruct.par.Option.Opt.Coord.CARTESIAN):
    lambda osp: 'set opt_coordinates cartesian',
    elstruct.option.name(elstruct.par.Option.Opt.Coord.ZMATRIX):
    lambda osp: 'set opt_coordinates internal',
    elstruct.option.name(elstruct.par.Option.Opt.Coord.REDUNDANT):
    lambda osp: 'set opt_coordinates redundant',
}
