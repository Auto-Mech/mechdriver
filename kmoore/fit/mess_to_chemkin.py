"""
  Routine for a single reaction:
    (1) Grab high-pressure and pressure rate constants
        from a MESS output file
    (2) Fit rate constants to an Arrhenius expression
    (3) Write Arrhenus parameters to a string formatted for
        a CHEMKIN mechanism file
"""

import arrfit
import mess_io
import chemkin_io


# MESS info
MESS_PATH = 'rate.out'
REACTANT = 'REACS'
PRODUCT = 'WR'
RXN = 'CH3+H=CH4'

# Desired info
FIT_PS = ['1', '10', '100', '1000', 'high']
TMIN = None
TMAX = None

# Dictionaries to store info; indexed by pressure (given in FIT_PS)
CALC_K_DCT = {}
FILT_CALC_TK_DCT = {}
FIT_PARAM_DCT = {}
FIT_K_DCT = {}
FIT_ERR_DCT = {}

# Fit info
T_REF = 298.0
FIT_TYPES = ['single', 'double']
FIT_METHOD = 'python'


# Loop over the single and double Arrhenius fits
for fit_type in FIT_TYPES:

    # Read the MESS output file into a string
    with open(MESS_PATH, 'r') as mess_file:
        output_string = mess_file.read()

    # Read the temperatures and pressures out of the MESS output
    mess_temps, _ = mess_io.reader.rates.get_temperatures(
        output_string)
    mess_pressures, punit = mess_io.reader.rates.get_pressures(
        output_string)

    # Loop over the pressures obtained from the MESS output
    for pressure in FIT_PS:

        assert pressure in mess_pressures

        # Read the rate constants
        if pressure == 'high':
            rate_ks = mess_io.reader.highp_ks(
                output_string, REACTANT, PRODUCT)
        else:
            rate_ks = mess_io.reader.pdep_ks(
                output_string, REACTANT, PRODUCT, pressure, punit)

        # Store in a the dictionary
        CALC_K_DCT[pressure] = rate_ks

    # Filter temperatures and rate_constants stored in the dictionary
    for pressure, calc_ks in CALC_K_DCT.items():
        flt_temps, flt_ks = arrfit.fit.get_valid_temps_rate_constants(
            mess_temps, calc_ks, tmin=TMIN, tmax=TMAX)

        # Store in a the dictionary
        FILT_CALC_TK_DCT[pressure] = [flt_temps, flt_ks]

    # Calculate the fitting parameters from the filtered T,k lists
    for pressure, tk_lsts in FILT_CALC_TK_DCT.items():

        # Set the temperatures and rate constants
        temps = tk_lsts[0]
        rate_constants = tk_lsts[1]

        # Obtain the fitting parameters based on the desired fit
        if fit_type == 'single' and FIT_METHOD == 'python':
            fit_params = arrfit.fit.single_arrhenius_fit(
                temps, rate_constants, T_REF)
        elif fit_type == 'double' and FIT_METHOD == 'python':
            init_params = arrfit.fit.single_arrhenius_fit(
                temps, rate_constants, T_REF)
            fit_params = arrfit.fit.double_arrhenius_fit_scipy(
                init_params[0], init_params[1], init_params[2],
                temps, rate_constants, T_REF)

        # Store the fitting parameters in a dictionary
        FIT_PARAM_DCT[pressure] = fit_params

    # Calculate fitted rate constants using the fitted parameters
    for pressure, params in FIT_PARAM_DCT.items():

        # Set the temperatures
        temps = FILT_CALC_TK_DCT[pressure][0]

        # Calculate fitted rate constants, based on fit type
        if fit_type == 'single':
            fit_ks = arrfit.fit.single_arrhenius(
                params[0], params[1], params[2],
                T_REF, temps)
        elif fit_type == 'double':
            fit_ks = arrfit.fit.double_arrhenius(
                params[0], params[1], params[2],
                params[3], params[4], params[5],
                T_REF, temps)

        # Store the fitting parameters in a dictionary
        FIT_K_DCT[pressure] = fit_ks

    # Calculate the error between the calc and fit ks
    for pressure, fit_ks in FIT_K_DCT.items():

        calc_ks = FILT_CALC_TK_DCT[pressure][1]
        sse, mean_avg_err, max_avg_err = arrfit.fit.calc_sse_and_mae(
            calc_ks, fit_ks)

        # Store in a dictionary
        FIT_ERR_DCT[pressure] = [sse, mean_avg_err, max_avg_err]

    # OUTPUT #

    # Print the fit type
    # print('\n\nFit type =   ', fit_type)

    # for pressure in FIT_PS:
    #     # Print the fit type
    #     print('\n\n  Pressure =   ', pressure)
    #     # Print the mess_tempsred t,k pairs type
    #     print('\n\n    Filtered T,K =   ', FILT_CALC_TK_DCT[pressure])
    #     # Print the fitting parameters
    #     print('\n    Fitting Parameters =   ', FIT_PARAM_DCT[pressure])
    #     # Print the fitting parameters
    #     print('\n    Rate Constants from fit =   ', FIT_K_DCT[pressure])
    #     # Print the fitting errors
    #     print('\n    Errors from fit =   ', FIT_ERR_DCT[pressure])

    # Write and print the pressure string
    pressure_str = chemkin_io.pfit.write_params(RXN, FIT_PARAM_DCT)
    print(pressure_str)
