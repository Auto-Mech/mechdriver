"""
  Library for writing output in CHEMKIN formatted files
"""

import os


def model_header(pf_levels, pf_models, refscheme = ''):
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
    if refscheme:    
        chemkin_header_str += '! reference scheme: {0}\n'.format(refscheme)

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


def write_rxn_file(ckin_rxn_dct, pes_formula, ckin_path):
    """ write out the rates
    """

    # Ensure the ckin directory path exists
    if not os.path.exists(ckin_path):
        os.makedirs(ckin_path)

    # Set the header string
    header_str = ckin_rxn_dct['header'] + '\n\n\n'

    # Write the header to the new, full PES ckin file
    # This also causes old file to be overwritten
    pes_ckin_name = os.path.join(ckin_path, pes_formula+'.ckin')
    with open(pes_ckin_name, 'w') as cfile:
        cfile.write(header_str)

    # Write the files by looping over the string
    for rxn, rstring in ckin_rxn_dct.items():
        if rxn != 'header':

            # Append to PES file
            with open(pes_ckin_name, 'a') as cfile:
                cfile.write(rstring+'\n\n')

            # Write to individual ckin file
            rxn_ckin_name = os.path.join(ckin_path, rxn+'.ckin')
            with open(rxn_ckin_name, 'w') as cfile:
                cfile.write(header_str + rstring)


def nasa_polynomial(hform0, hform298, ckin_poly_str):
    """ write the nasa polynomial str
    """
    hf_str = (
        '! Hf(0 K) = {:.2f},'.format(hform0) +
        'Hf(298 K) = {:.2f} kcal/mol\n'.format(hform298)
    )
    return hf_str + ckin_poly_str


def write_nasa_file(ckin_nasa_str, ckin_path):
    """ write out the nasa polynomials
    """
    if not os.path.exists(ckin_path):
        os.makedirs(ckin_path)
    fpath = os.path.join(ckin_path, 'all.ckin')
    with open(fpath, 'w') as nasa_file:
        nasa_file.write(ckin_nasa_str)


def write_transport_file(ckin_trans_str, ckin_path):
    """ write out the transport file
    """
    if not os.path.exists(ckin_path):
        os.makedirs(ckin_path)
    fpath = os.path.join(ckin_path, 'trans.ckin')
    with open(fpath, 'w') as nasa_file:
        nasa_file.write(ckin_trans_str)
