""" test arrfit
"""

import chemkin_io
import scripts

PAIRS = [
    [600, 93.5678],
    [700, 2287.02],
    [800, 27103.4],
    [900, 193348],
    [1000, 955781],
    [1100, 3.60E+06],
    [1200, 1.10E+07],
    [1300, 2.85E+07],
    [1400, 6.51E+07],
    [1500, 1.34E+08],
    [1600, 2.52E+08],
    [1700, 4.43E+08],
    [1800, 7.34E+08],
    [1900, 1.16E+09],
    [2000, 1.74E+09],
    [2100, 2.53E+09],
    [2200, 3.56E+09],
    [2300, 4.86E+09],
    [2400, 6.48E+09],
    [2500, 8.45E+09],
    [2600, 1.08E+10],
    [2700, 1.36E+10],
    [2800, 1.68E+10],
    [2900, 2.06E+10],
    [3000, 2.48E+10]
]

MESS_PATH = '.'
REACTION = 'CH3+H=CH4'
ERR_THRESH = 2

KTP_DCT = {
#     1:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
#     10:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
#     100:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
    'high':  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
}


# Fit rate constants to single Arrhenius expressions
sing_rate_params, sing_errs = scripts.ktp.mod_arr_fit(
    KTP_DCT, MESS_PATH, fit_type='single', fit_method='dsarrfit',
    t_ref=1.0, a_conv_factor=1.0)

# Assess error from single fitting and ntemps to decide if dbl fit to be done
sgl_fit_good = max((vals[1] for vals in sing_errs.values())) < ERR_THRESH
dbl_fit_imposs = any(len(KTP_DCT[p][0]) < 6 for p in KTP_DCT)

# Write chemkin string for single, or perform dbl fit and write string
chemkin_str = ''
if sgl_fit_good:
    if not dbl_fit_imposs:
        print('Single fit errors acceptable: Using single fits')
    if dbl_fit_imposs:
        print('Not enough temperatures for a double fit: Using single fits')
    chemkin_str += chemkin_io.writer.reaction.plog(
        REACTION, sing_rate_params, sing_errs)
else:
    print('Single fit errors too large and double fitting possible: Doing a double fit')
    doub_rate_params, doub_errs = scripts.ktp.mod_arr_fit(
        KTP_DCT, MESS_PATH, fit_type='double', fit_method='dsarrfit',
        t_ref=1.0, a_conv_factor=1.0)
    chemkin_str += chemkin_io.writer.reaction.plog(
        REACTION, doub_rate_params, doub_errs)

print(chemkin_str)
