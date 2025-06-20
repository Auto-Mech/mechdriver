""" test mechanalyzer pssa tools 
    combined builder / calculator bf functions
"""

import os
import copy
import numpy as np
from ioformat import pathtools, remove_comment_lines
from mechanalyzer.parser import ckin_ as ckin_parser
from mechanalyzer.parser._util import resort_ktp_labels
from mechanalyzer import calculator
import chemkin_io # better to remove, but would need too many data to write in this file

CWD = os.path.dirname(os.path.realpath(__file__))
FILE1 = 'data/pssa/C6H5C3H3.mech'
FILE2 = 'data/pssa/C7H5C3H4.mech'
THERM1 = 'data/pssa/C6H5C3H3.therm'
PSSA_SPCS1 = ['C6H5C3H3']
PSSA_SPCS2 = ['BENZYLCPDYL-M+H','p6+H','p4+H','p16+H','p10+H','p13+C2H2', 'mindrC','meindrB','AZLW','FLW','FULVALENE+H','AZULENE+H' ]
TVECT = [(np.arange(500, 2100, 100))]
PVECT = [0.1, 1., 10., 100]

# PARTIAL RESULTS, SELECTED AT RANDOM P AND RANDOM RXNS
RESULTS1 = {(('INDENYL', 'H'), ('INDENE',), (None,)): {0.1: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([1.63930489e+14, 2.72530235e+14, 3.34443865e+14, 3.46328923e+14,
       3.24525304e+14, 2.86184152e+14, 2.43102846e+14, 2.01805563e+14,
       1.65220786e+14, 1.34212675e+14, 1.08606608e+14, 8.77855333e+13,
       7.10047340e+13, 5.75424043e+13, 4.67614865e+13, 3.81267001e+13]))}, 
       (('INDENYL', 'H'), ('C6H5', 'C3H3'), (None,)): {100: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([2.57818544e-24, 1.48296172e-17, 4.48245914e-13, 7.25208817e-10,
       2.40432453e-07, 2.81511062e-05, 1.55196422e-03, 4.83105292e-02,
       9.61697180e-01, 1.33838930e+01, 1.39067968e+02, 1.13373067e+03,
       7.53282128e+03, 4.20157413e+04, 2.01340418e+05, 8.44148938e+05]))}, 
       (('INDENE',), ('INDENYL', 'H'), (None,)): {1.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([4.42615882e-21, 8.65283173e-15, 2.33164564e-10, 4.39779924e-07,
       1.42609379e-04, 1.35987060e-02, 5.35723889e-01, 1.09261196e+01,
       1.34771721e+02, 1.12292451e+03, 6.85120041e+03, 3.25094089e+04,
       1.25600298e+05, 4.09380125e+05, 1.15752425e+06, 2.90298806e+06]))}, 
       (('C6H5', 'C3H3'), ('INDENYL', 'H'), (None,)): {10.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([1.20463990e+03, 5.16126421e+04, 9.44542140e+05, 9.87636564e+06,
       6.97925956e+07, 3.70037474e+08, 1.57697397e+09, 5.66473682e+09,
       1.77448410e+10, 4.97032058e+10, 1.26866184e+11, 2.99443432e+11,
       6.61143438e+11, 1.37813473e+12, 2.73243865e+12, 5.18489419e+12]))}, 
       (('INDENE',), ('C6H5', 'C3H3'), (None,)): {1.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([7.08135949e-31, 2.32014390e-22, 1.80862472e-16, 3.29810416e-12,
       5.09900551e-09, 1.46255808e-06, 1.27578442e-04, 4.67224685e-03,
       8.92417094e-02, 1.03394640e+00, 8.10937072e+00, 4.68848271e+01,
       2.15288558e+02, 8.48500656e+02, 3.13357699e+03, 1.17336203e+04])), 10.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([5.09880660e-30, 1.54105145e-21, 8.58049051e-16, 1.01012254e-11,
       1.13521691e-08, 2.75427401e-06, 2.27496819e-04, 8.49852356e-03,
       1.73620352e-01, 2.21508989e+00, 1.94321298e+01, 1.26084849e+02,
       6.42801272e+02, 2.73840866e+03, 1.05913414e+04, 4.17778807e+04]))}}
RESULTS2 = {(('C7H5', 'C3H4-A'), ('C9H6CH2', 'H'), (None,)): {0.1: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([2.00852107e+04, 2.04319512e+05, 1.04196758e+06, 3.46326345e+06,
       8.67328220e+06, 1.78456630e+07, 3.18655792e+07, 5.12064054e+07,
       7.59278247e+07, 1.05746824e+08, 1.40136260e+08, 1.78421697e+08,
       2.19862227e+08, 2.63710958e+08, 3.09256261e+08, 3.55847105e+08]))}, 
       (('MEINDENYL',), ('C7H5', 'C3H4-A'), (None,)): {1.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([7.83297579e-45, 5.17644208e-34, 1.99020169e-26, 7.58307794e-21,
       1.37051366e-16, 2.98944594e-13, 1.42077977e-10, 2.17624993e-08,
       1.40648301e-06, 4.64376090e-05, 9.00452274e-04, 1.13773153e-02,
       1.01369898e-01, 6.76991570e-01, 3.55521349e+00, 1.52508867e+01]))},
       (('C7H5', 'C3H4-A'), ('NAPH', 'H'), (None,)): {10.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([9.39301208e+04, 4.36291799e+05, 1.33005767e+06, 3.13429994e+06,
       6.33070998e+06, 1.14098892e+07, 1.89383194e+07, 3.10213968e+07,
       5.16346133e+07, 8.42968712e+07, 1.30300322e+08, 1.89084647e+08,
       2.58464986e+08, 3.33995035e+08, 4.09423245e+08, 4.78399725e+08]))}, 
       (('C7H5', 'C3H4-A'), ('C7H5', 'C3H4-P'), (None,)): {0.1: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([4.97139190e-05, 3.15620241e-02, 4.61914397e+00, 2.07466404e+02,
       3.58709427e+03, 3.20882427e+04, 1.88270556e+05, 8.32565484e+05,
       3.01230371e+06, 9.42563585e+06, 2.65319026e+07, 6.86741631e+07,
       1.63812860e+08, 3.57260794e+08, 7.10886511e+08, 1.30221734e+09]))}, 
       (('MEINDENYL',), ('C10H9',), (None,)): {100: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([4.16423607e-31, 3.76981796e-22, 4.06700521e-16, 7.25720798e-12,
       9.01890718e-09, 1.82342520e-06, 1.01987455e-04, 2.23487736e-03,
       2.43131086e-02, 1.55052931e-01, 6.53374879e-01, 1.98691263e+00,
       4.65879915e+00, 8.85907570e+00, 1.42073926e+01, 1.98144802e+01]))},
       (('MEINDRSR',), ('C10H9',), (None,)): {1.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([5.55951945e-15, 3.10739061e-08, 5.92599427e-04, 3.79452134e-01,
       2.80565938e+01, 4.92239287e+02, 3.19785794e+03, 1.02589340e+04,
       1.97154923e+04, 2.59427335e+04, 2.56983472e+04, 2.05268636e+04,
       1.39091007e+04, 8.30526724e+03, 4.49822670e+03, 2.25952627e+03]))}, 
       (('C9H6CH2', 'H'), ('MEINDENYL',), (None,)): {10.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([5.87108819e-02, 2.80849243e+02, 8.40933302e+04, 4.60716670e+06,
       8.37545057e+07, 7.19144879e+08, 3.63971947e+09, 1.25599175e+10,
       3.26401970e+10, 6.85229681e+10, 1.22278259e+11, 1.92629974e+11,
       2.75680760e+11, 3.66399516e+11, 4.60005077e+11, 5.52776463e+11]))}, 
       (('NAPH', 'H'), ('C9H6CH2', 'H'), (None,)): {1.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([1.75770368e-03, 1.73764935e+01, 9.11897307e+03, 7.94070612e+05,
       2.14243372e+07, 2.59156339e+08, 1.77219734e+09, 7.97874647e+09,
       2.62414417e+10, 6.78357612e+10, 1.45316708e+11, 2.68242186e+11,
       4.39413127e+11, 6.53341877e+11, 8.97295826e+11, 1.15410959e+12]))}, 
       (('C7H5', 'C3H4-P'), ('MEINDENYL',), (None,)): {0.1: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([1.13402193e-03, 1.16655264e+00, 1.06232967e+02, 2.26558008e+03,
       1.88831105e+04, 8.25943313e+04, 2.28582855e+05, 4.54476657e+05,
       7.09694449e+05, 9.27464091e+05, 1.06162278e+06, 1.09992735e+06,
       1.05645129e+06, 9.57270470e+05, 8.29029021e+05, 6.92965894e+05]))}, 
       (('C10H9',), ('NAPH', 'H'), (None,)): {10.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([2.57264533e-06, 1.49054550e-01, 1.66380656e+02, 1.74374093e+04,
       4.04205077e+05, 3.41790223e+06, 1.43696842e+07, 3.67118393e+07,
       6.52265020e+07, 8.84926996e+07, 9.79770620e+07, 9.29061264e+07,
       7.81954779e+07, 6.00064568e+07, 4.28529183e+07, 2.89322972e+07])), 100: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([2.57264533e-06, 1.49054550e-01, 1.66380656e+02, 1.74374093e+04,
       4.04205077e+05, 3.41790223e+06, 1.43696842e+07, 3.67118393e+07,
       6.52265020e+07, 8.84926996e+07, 9.79770620e+07, 9.29061264e+07,
       7.81954779e+07, 6.00064568e+07, 4.28529183e+07, 2.89322972e+07]))},
       (('MEINDRSR',), ('C9H6CH2', 'H'), (None,)): {0.1: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([1.16293523e-07, 2.22087942e-03, 1.40123624e+00, 1.13072114e+02,
       2.43418303e+03, 2.15146842e+04, 1.02068259e+05, 3.09442841e+05,
       6.74507055e+05, 1.14752979e+06, 1.61583833e+06, 1.96564907e+06,
       2.13278037e+06, 2.11447249e+06, 1.95127163e+06, 1.70036616e+06])), 1.0: (np.array([ 500,  600,  700,  800,  900, 1000, 1100, 1200, 1300, 1400, 1500,
       1600, 1700, 1800, 1900, 2000]), np.array([1.38921028e-08, 6.76258021e-04, 8.24341696e-01, 1.08209951e+02,
       3.38175462e+03, 4.00906438e+04, 2.40941276e+05, 8.86837057e+05,
       2.27190202e+06, 4.42901856e+06, 7.00402750e+06, 9.41497250e+06,
       1.11395009e+07, 1.19115267e+07, 1.17476378e+07, 1.08564330e+07]))}}

#####################################################################

# ktp dct

def test_pssa_wtherm():
    rxn_ktp_dct = read_inputs(FILE1)
    ######################################## COMPUTE BW RXNS
    rxn_ktp_dct_bw = {}
    # read thermo
    therm_str = pathtools.read_file(CWD, THERM1, remove_comments='!')
    spc_therm_dct = ckin_parser.parse_spc_therm_dct(therm_str, TVECT[0])
    spc_term_df = calculator.thermo.spc_therm_dct_df(spc_therm_dct)
    for rxn, k in rxn_ktp_dct.items():
        newlabel = (rxn[1], rxn[0], rxn[2])
        if newlabel not in rxn_ktp_dct.keys():
            newentry = calculator.rates.get_bw_ktp_dct(k, rxn[0], rxn[1], spc_term_df)
            rxn_ktp_dct_bw[newlabel] = newentry
    # merge
    rxn_ktp_dct = calculator.rates.merge_rxn_ktp_dcts(
        rxn_ktp_dct, rxn_ktp_dct_bw)
    ############################################################
    # GET KTP DICTIONARY
    rxn_ktp_dct = get_pssa_dct(rxn_ktp_dct, PSSA_SPCS1)

    # assess results
    for rxn, rxn_dct in RESULTS1.items():
        for p, temp_and_k in rxn_dct.items():
            arr1 = temp_and_k[1]
            arr2 = rxn_ktp_dct[rxn][p][1]
            assert np.allclose(arr1, arr2, atol=1e-3)

def test_pssa_multiplespc():
    rxn_ktp_dct = read_inputs(FILE2)
    rxn_ktp_dct = get_pssa_dct(rxn_ktp_dct, PSSA_SPCS2)

    # assess results
    for rxn, rxn_dct in RESULTS2.items():
        for p, temp_and_k in rxn_dct.items():
            arr1 = temp_and_k[1]
            arr2 = rxn_ktp_dct[rxn][p][1]
            assert np.allclose(arr1, arr2, atol=1e-3)


def read_inputs(inpfile):
    # read mech AND ORDER NAMES OF RCTS AND PRDS IN THE SAME WAY
    cki_out = pathtools.read_file(CWD, inpfile)
    cki_out = remove_comment_lines(cki_out, '!')
    rxn_param_dct = chemkin_io.parser.reaction.get_rxn_param_dct(
        cki_out, 'cal/mole', 'moles')
    rxn_param_dct = resort_ktp_labels(rxn_param_dct)
    rxn_ktp_dct = calculator.rates.eval_rxn_param_dct(
        rxn_param_dct, TVECT, PVECT)
    return rxn_ktp_dct

def get_pssa_dct(rxn_ktp_dct, pssa_spc):
    for SP_TOREPLACE in pssa_spc:
        rxn_ktp_dct_new = {} # new ktp dct
        # get bfs for species
        bf_df_sp = calculator.bf.bf_df_fromktpdct(rxn_ktp_dct, SP_TOREPLACE, TVECT[0], PVECT)
        bf_dct = calculator.bf.bf_tp_df_todct(bf_df_sp, bf_threshold = 0.1)
        for rxn, k in rxn_ktp_dct.items():
            if tuple(SP_TOREPLACE.split('+')) == rxn[1]:
                # species in products: replace corresponding branching fractions
                for prd, bf_sp in bf_dct.items():
                    # write a new reaction having the new species as product
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

    return rxn_ktp_dct

if __name__ == '__main__':
    test_pssa_wtherm()
    test_pssa_multiplespc()
