
"""
Apply PSSA on one or more species on a PES
works with input in mechanalyzer/tests/data/pssa/C7H5C3H4.mech
"""

import os
import copy
import numpy as np
from ioformat import pathtools, remove_comment_lines
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.parser._util import resort_ktp_labels
from mechanalyzer import calculator
import chemkin_io
import ratefit


def main(
    startmech: str='kin.CKI',
    pssa_spcs: list=[],
    temprange: list = [500, 2000],
    pvect: list = [1.0],
    bfthresh: float = 1e-4,
    thermofile: str = None,
    outputrates: str = 'pssa_rates.txt',
):
    Tvect = [np.arange(temprange[0], temprange[1]+100, 100)]
    pvect = list(pvect)
    #####################################################################
    # Set path to current directory where MESS files exist
    CWD = os.getcwd()

    # read mech AND ORDER NAMES OF RCTS AND PRDS IN THE SAME WAY
    cki_out = pathtools.read_file(CWD, startmech)
    cki_out = remove_comment_lines(cki_out, '!')
    rxn_param_dct = chemkin_io.parser.reaction.get_rxn_param_dct(
        cki_out, 'cal/mole', 'moles')
    rxn_param_dct = resort_ktp_labels(rxn_param_dct)
    rxn_ktp_dct = calculator.rates.eval_rxn_param_dct(
        rxn_param_dct, Tvect, pvect)

    # HERE: COMPUTE THE RATE CONSTANT FOR THE BACKWARD REACTION
    if thermofile:
        rxn_ktp_dct_bw = {}
        # read thermo
        therm_str = pathtools.read_file(CWD, thermofile, remove_comments='!')
        spc_therm_dct = ckin_parser.parse_spc_therm_dct(therm_str, Tvect[0])
        spc_term_df = calculator.thermo.spc_therm_dct_df(spc_therm_dct)
        for rxn, k in rxn_ktp_dct.items():
            newlabel = (rxn[1], rxn[0], rxn[2])
            if newlabel not in rxn_ktp_dct.keys():
                newentry = calculator.rates.get_bw_ktp_dct(k, rxn[0], rxn[1], spc_term_df)
                rxn_ktp_dct_bw[newlabel] = newentry
        # merge
        rxn_ktp_dct = calculator.rates.merge_rxn_ktp_dcts(
            rxn_ktp_dct, rxn_ktp_dct_bw)
    
    rxn_ktp_dct_new = {}
    for SP_TOREPLACE in pssa_spcs:
        rxn_ktp_dct_new = {} # new ktp dct
        # get bfs for species
        bf_df_sp = calculator.bf.bf_df_fromktpdct(rxn_ktp_dct, SP_TOREPLACE, Tvect[0], pvect)
        bf_dct = calculator.bf.bf_tp_df_todct(bf_df_sp, bf_threshold = bfthresh)
        for rxn, k in rxn_ktp_dct.items():
            if tuple(SP_TOREPLACE.split('+')) == rxn[1]:
                # species in products: replace corresponding branching fractions
                for prd, bf_sp in bf_dct.items():
                    # write a new reaction having the new species as product
                    print({rxn: k}, bf_sp)
                    ktp_sp = calculator.bf.merge_bf_rates({rxn: k}, bf_sp)
                    newlabel = (rxn[0], tuple(prd.split('+')), rxn[2])   
                    rxn_ktp_dct_new = calculator.rates.merge_rxn_ktp_dcts(
                        rxn_ktp_dct_new, {newlabel: ktp_sp[rxn]})
    
            elif tuple(SP_TOREPLACE.split('+')) == rxn[0]:
                continue # the species is a reactant - do not add the reaction
            else:
                # add the reaction with the rest
                rxn_ktp_dct_new = calculator.rates.merge_rxn_ktp_dcts(
                    rxn_ktp_dct_new, {rxn: k})  
        # the new ktp dct will become the rxt_ktp_dct from which you will extract new product branching fractions
        rxn_ktp_dct = copy.deepcopy(rxn_ktp_dct_new)

    # Generate CKIN files
    # Fit
    rxn_param_dct, rxn_err_dct = ratefit.fit.fit_rxn_ktp_dct(
        rxn_ktp_dct_new, 'arr', arrfit_dct={'dbltol': 1.},
        pdep_dct={'temps': (500, 1000, 2000), 'tol': 0.01,
                'plow': None, 'phigh': None, 'pval': 1.0}
    )

    # Get the comments dct and write the Chemkin string
    rxn_cmts_dct = chemkin_io.writer.comments.get_rxn_cmts_dct(
        rxn_err_dct=rxn_err_dct)
    ckin_str = chemkin_io.writer.mechanism.write_chemkin_file(
        rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct, sortrxns=True)

    # Write the fitted rate parameter file
    pathtools.write_file(ckin_str, CWD, outputrates)

