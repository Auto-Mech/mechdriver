""" test arrfit
"""
import numpy
import chemkin_io
import scripts


PAIRS2 = [
    [700, 2287.02],
    [900, 193348],
    [1000, 955781],
    [1100, 3.60E+06]
]
PAIRS3 = [
    [1900, 1.16E+09],
    [2000, 1.74E+09],
    [2100, 2.53E+09],
    [2200, 3.56E+09],
    [2300, 4.86E+09],
]
PAIRS = [
#    [600, 93.5678],
#    [700, 2287.02],
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
ERR_THRESH = 200

KTP_DCT = {
#     1:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
#     10:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
#      10:  [[pair[0] for pair in PAIRS2], [pair[1] for pair in PAIRS2]],
      10:  numpy.array([[pair[0] for pair in PAIRS2], [pair[1] for pair in PAIRS2]]),
#     100:  [[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]],
     'high':  numpy.array([[pair[0] for pair in PAIRS], [pair[1] for pair in PAIRS]]),
}


# Fit rate constants to single Arrhenius expressions
sing_params_dct, sing_fit_success = scripts.ktp.mod_arr_fit(
    KTP_DCT, MESS_PATH, fit_type='single', fit_method='dsarrfit',
    t_ref=1.0, a_conv_factor=1.0)
if sing_fit_success:
    print('Successful fit to Single Arrhenius at all T, P')
    print(sing_params_dct)

# Assess the errors of the single Arrhenius Fit
sing_fit_err_dct = scripts.ktp.assess_arr_fit_err(
    sing_params_dct, KTP_DCT, fit_type='single',
    t_ref=1.0, a_conv_factor=1.0)

# Assess error from single fitting and ntemps to decide if dbl fit to be done
sgl_fit_good = max((vals[1] for vals in sing_fit_err_dct.values())) < ERR_THRESH
dbl_fit_poss = any(len(KTP_DCT[p][0]) >= 6 for p in KTP_DCT)

# Write chemkin string for single, or perform dbl fit and write string
chemkin_str = ''
if sgl_fit_good:
    print('Single fit errors acceptable: Using single fits')
    chemkin_str += chemkin_io.writer.reaction.plog(
        REACTION, sing_params_dct, sing_fit_err_dct)
elif not sgl_fit_good and dbl_fit_poss:
    print('Single fit errs too large & double fit possible: Trying double fit')

    # Fit rate constants to double Arrhenius expressions
    doub_params_dct, doub_fit_success = scripts.ktp.mod_arr_fit(
        KTP_DCT, MESS_PATH, fit_type='double', fit_method='dsarrfit',
        t_ref=1.0, a_conv_factor=1.0)

    if doub_fit_success:
        print('Successful fit to Double Arrhenius at all T, P')
        # Assess the errors of the single Arrhenius Fit
        doub_fit_err_dct = scripts.ktp.assess_arr_fit_err(
            doub_params_dct, KTP_DCT, fit_type='double',
            t_ref=1.0, a_conv_factor=1.0)
        chemkin_str += chemkin_io.writer.reaction.plog(
           REACTION, doub_params_dct, doub_fit_err_dct)
    else:
        print('Double Arrhenius Fit failed for some reason: Using Single fits')
        chemkin_str += chemkin_io.writer.reaction.plog(
            REACTION, sing_params_dct, sing_fit_err_dct)
elif not sgl_fit_good and not dbl_fit_poss:
        print('Not enough temperatures for a double fit: Using single fits')
        chemkin_str += chemkin_io.writer.reaction.plog(
            REACTION, sing_params_dct, sing_fit_err_dct)

print(chemkin_str)
