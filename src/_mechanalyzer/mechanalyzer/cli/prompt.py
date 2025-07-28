
"""
Analyze the extent of prompt dissociation
for a single exothermic reaction with successive decomposition
"""

import os
from ioformat import pathtools, remove_comment_lines
import autoparse.pattern as app
from mess_io import reader
from chemkin_io import writer
from mechanalyzer.builder import multipes_prompt_dissociation_ktp_dct
import ratefit

CWD = os.getcwd()
def main(
    flds: list = ['',],
    messinput: str = 'mess.inp',
    messoutput: str = 'mess.out',
    outputmicro: str = 'ke.out',
    log: str = 'mess.log',
    model: str = 'rovib_dos',
    bfthresh: float = 0.1,
    fitmethod: str = 'plog',
    outputrates: str = 'rates_prompt.txt'
):

    list_strs_dct = []
    # Read the input and output files for MESS calculation of 1st PES
    # if no directories are listed: add all directories from the current one
    if (flds[0] == '' and len(flds) == 1) or len(flds) == 0:
        print('*Warning: no directories given- all dirs in current dir will be considered')
        flds = [d for d in os.listdir(CWD) if os.path.isdir(os.path.join(d))]
    # Here Product Energy Distributions are calculated
    for FLD in flds:
        flddct = dict.fromkeys(['inp', 'ktp_out', 'ke_out', 'ped', 'log'])
        FLDPATH = os.path.join(CWD, FLD)
        me_inp = pathtools.read_file(FLDPATH, messinput)
        me_inp = remove_comment_lines(me_inp, delim_pattern=app.escape('!'))
        me_inp = remove_comment_lines(me_inp, delim_pattern=app.escape('#'))
        flddct['inp'] = me_inp
        flddct['ktp_out'] = pathtools.read_file(FLDPATH, messoutput)
        flddct['ke_out'] = pathtools.read_file(FLDPATH, outputmicro)
        _, pedoutput = reader.ped.ped_names(me_inp)
        if pedoutput:
            flddct['ped'] = pathtools.read_file(FLDPATH, pedoutput)
        flddct['log'] = pathtools.read_file(FLDPATH, log)
        list_strs_dct.append(flddct)

    # Read and calculate the rate constants into a rxn ktp dictionary
    rxn_ktp_dct = multipes_prompt_dissociation_ktp_dct(
        list_strs_dct,
        model, bfthresh
    )

    # Fit
    rxn_param_dct, rxn_err_dct = ratefit.fit.fit_rxn_ktp_dct(
        rxn_ktp_dct, fitmethod, arrfit_dct={'dbltol': 50.},
        pdep_dct={'temps': (500, 1000, 2000), 'tol': 0.01,
                'plow': None, 'phigh': None, 'pval': 1.0}
    )

    # Get the comments dct and write the Chemkin string
    rxn_cmts_dct = writer.comments.get_rxn_cmts_dct(
        rxn_err_dct=rxn_err_dct)
    ckin_str = writer.mechanism.write_chemkin_file(
        rxn_param_dct=rxn_param_dct, rxn_cmts_dct=rxn_cmts_dct)

    # Write the fitted rate parameter file
    pathtools.write_file(ckin_str, CWD, outputrates)
