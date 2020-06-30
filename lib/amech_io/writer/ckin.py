"""
  Library for writing output in CHEMKIN formatted files
"""

import os
import automol
import routines.pf.thermo


def model_header(pf_levels, pf_models):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """

    tors_model = pf_models['tors']
    vib_model = pf_models['vib']
    sym_model = pf_models['sym']
    geo_info = pf_levels['geo'][1]
    ene_info = pf_levels['ene'][1]
    har_info = pf_levels['harm'][1]
    tors_info = pf_levels['tors'][1]

    # Convert the pac99 polynomial to chemkin polynomial
    chemkin_header_str = '! vib model: {0}\n'.format(vib_model)
    chemkin_header_str += '! tors model: {0}\n'.format(tors_model)
    # chemkin_header_str += '! vpt2 model: {0}\n'.format(vpt2_model)
    chemkin_header_str += '! sym model: {0}\n'.format(sym_model)
    if har_info and ene_info:
        chemkin_header_str += _ckin_ene_lvl_str(ene_info, geo_info)
    if tors_info:
        chemkin_header_str += '! tors level: {}{}/{}//{}{}/{}\n'.format(
            tors_info[1][3], tors_info[1][1], tors_info[1][2],
            tors_info[0][3], tors_info[0][1], tors_info[0][2])
    # if vpt2_info:
    #     chemkin_header_str += '! vpt2 level: {}{}/{}\n'.format(
    #         vpt2_info[3], vpt2_info[1], vpt2_info[2])

    return chemkin_header_str


def _ckin_ene_lvl_str(ene_info, geo_info):
    """ Write the comment lines for the enrgy lvls for ckin
    """
    ene_str = '! energy level:'
    for i, ene_lvl in enumerate(ene_info):
        ene_str += ' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
            ene_lvl[0],
            ene_lvl[1][3], ene_lvl[1][1], ene_lvl[1][2],
            geo_info[3], geo_info[1], geo_info[2])
        if i+1 != len(ene_info):
            ene_str += ' +'

    return ene_str


def run_ckin_poly(spc, spc_dct_i, pac99_poly_str):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """

    hf_str = '! Hf(0 K) = {:.2f}, Hf(298 K) = {:.2f} kcal/mol\n'.format(
        float(spc_dct_i['Hfs'][0]), float(spc_dct_i['Hfs'][1]))
    ich = spc_dct_i['ich']
    formula_dct = automol.inchi.formula(ich)
    chemkin_poly_str = routines.pf.thermo.nasapoly.convert_pac_to_chemkin(
        spc, formula_dct, hf_str, pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)
    return chemkin_poly_str


def write_rates_file(ckin_rate_str_lst):
    """ write out the rates
    """

    # Write header string containing thy information
    chemkin_header_str = writer.ckin.model_header(es_info, pf_model)
    chemkin_header_str += '\n'

    # Initialize full chemkin string and paths
    # chemkin_full_str = chemkin_header_str
    chemkin_full_str = ''
    starting_path = os.getcwd()
    ckin_path = ''.join([starting_path, '/ckin'])
    if not os.path.exists(ckin_path):
        os.mkdir(ckin_path)
    full_ckin_path = os.path.join(ckin_path, pes_formula_str+'.ckin')
    # if os.path.exists(full_ckin_path):
    #     os.remove(full_ckin_path)


    ####
                    # Print the results for each channel to a file
                    chemkin_full_str += '\n\n'
                    pes_chn_lab = pes_formula_str + '_' + name_i + '_' + name_j
                    ckin_name = os.path.join(ckin_path, pes_chn_lab+'.ckin')
                    with open(ckin_name, 'w') as cfile:
                        cfile.write(chemkin_str)

    # Print the results for the whole PES to a file
    with open(full_ckin_path, 'a') as cfile:
        cfile.write(chemkin_full_str)


def write_nasa_file(ckin_path, chemkin_poly_str):
    """ write out the nasa polynomials
    """
    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
