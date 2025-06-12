""" elstruct.writer._molpro2015 parameters
"""
import sys
from elstruct import Option
from elstruct import option
import elstruct.par
from elstruct.writer import fill


REF_DCT = {
    elstruct.par.Reference.RHF: 'rhf',
    elstruct.par.Reference.UHF: 'uhf'
}


class MultiReference():
    """ _ """
    CASSCF = 'casscf'
    CASPT2 = 'rs2'
    CASPT2I = 'rs2'
    CASPT2C = 'rs2c'
    MRCI_Q = 'mrci'


OPTION_EVAL_DCT = {
    option.name(elstruct.par.Option.Mol.NOSYMM_):
    'nosym',
    option.name(Option.Scf.MAXITER_):
    lambda osp: 'maxit={}'.format(*option.values(osp)),
    option.name(Option.Scf.DIIS_):
    lambda osp: ('iptyp=diis' if option.values(osp)[0] else 'iptyp=none'),
    option.name(Option.Casscf.OCC_):
    lambda osp: 'occ,{}'.format(*option.values(osp)),
    option.name(Option.Casscf.CLOSED_):
    lambda osp: 'closed,{}'.format(*option.values(osp)),
    option.name(Option.Casscf.WFN_):
    lambda osp: 'wf,{},{},{},{};state,{}'.format(*option.values(osp)),
    option.name(Option.MRCorr.SHIFT_):
    lambda osp: 'shift={}'.format(*option.values(osp)),
    option.name(Option.Opt.MAXITER_):
    lambda osp: 'maxit={}'.format(*option.values(osp)),
    option.name(Option.Corr.ALL_ELEC_):
    lambda osp: '\n core,0',
}


def set_method_and_options(method, corr_options, gen_lines):
    """ Set corr options based on the method
    """

    core_method, pfxs = elstruct.par.Method.evaluate_method_type(method)
    
    # Add to options lists
    if elstruct.par.Method.ModPrefix.ALL_ELEC[0] in pfxs:
        allelec_opt = Option.Corr.ALL_ELEC_
        if allelec_opt not in corr_options:
            corr_options += (allelec_opt,)
    # gdirect here maybe as gen_lines? ONLY FOR DF-MP2 CASES
    if elstruct.par.Method.ModPrefix.REL_DKH[0] in pfxs:
        gen_lines = fill.update_gen_lines(
            gen_lines, lines2=('dkroll=1',))
        # corr_options, scf_options = ['cwrelativistic']
        # -> 'expec,rel,darwin
        # {method;
        #  expec,rel,darwin,massv}

    # Update the method that should be written to file
    fin_method = ''
    if elstruct.par.Method.ModPrefix.L_PNO[0] in pfxs:
        fin_method = method

    if elstruct.par.Method.ModPrefix.DF[0] in pfxs:
        fin_method = method
        gen_lines = fill.update_gen_lines(
            gen_lines, lines2=('gdirect',))       
         
    return core_method, fin_method, corr_options, gen_lines
