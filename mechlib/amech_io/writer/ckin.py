"""
  Library for writing output in CHEMKIN formatted files
"""

import os


# COMMON HEADER STUFF FOR kTP and THERMO CKIN FILES
def model_header(spc_mods, spc_mod_dct, refscheme=''):
    """ Write a model header for multiple models
    """
    mod_str = ''
    for spc_mod in spc_mods:
        mod_str += _model_header(spc_mod_dct[spc_mod], refscheme=refscheme)

    return mod_str


def _model_header(spc_mod_dct_i, refscheme=''):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """

    # Pull information out of pf dictionaries
    tors_model = spc_mod_dct_i['tors']['mod']
    vib_model = spc_mod_dct_i['vib']['mod']
    sym_model = spc_mod_dct_i['symm']['mod']

    har_info = spc_mod_dct_i['vib']['geolvl'][1]
    tors_geo_info = spc_mod_dct_i['tors']['geolvl'][1]
    tors_ene_info = spc_mod_dct_i['tors']['enelvl'][1]
    # vpt2_info = spc_mod_dct_i['vib']['vpt2lvl'][1]
    # vpt2_info = None

    ene_infos = tuple(inf for key, inf in spc_mod_dct_i['ene'].items()
                      if 'lvl' in key)

    # Write the theory information info header for a reaction/species
    chemkin_header_str = '! vib model: {0}\n'.format(vib_model)
    chemkin_header_str += '! tors model: {0}\n'.format(tors_model)
    chemkin_header_str += '! sym model: {0}\n'.format(sym_model)
    if har_info is not None and ene_infos:
        chemkin_header_str += _ckin_ene_lvl_str(ene_infos, har_info)
    if tors_geo_info is not None and tors_ene_info is not None:
        chemkin_header_str += '! tors level: {}{}/{}//{}{}/{}\n'.format(
            tors_ene_info[1][3], tors_ene_info[1][1], tors_ene_info[1][2],
            tors_geo_info[1][3], tors_geo_info[1][1], tors_geo_info[1][2])
    # if vpt2_info is not None:
    #     chemkin_header_str += '! vpt2 level: {}/{}\n'.format(
    #         vpt2_info[1][1], vpt2_info[1][2])
    if refscheme:
        chemkin_header_str += '! reference scheme: {0}\n'.format(refscheme)

    return chemkin_header_str


def _ckin_ene_lvl_str(ene_infos, har_info):
    """ Write the comment lines for the enrgy lvls for ckin
    """
    ene_str = '! energy level:'
    for i, ene_inf in enumerate(ene_infos):
        # print('ene inf test', ene_inf)
        ene_str += ' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
            ene_inf[1][0],
            ene_inf[1][1][3], ene_inf[1][1][1], ene_inf[1][1][2],
            har_info[1][3], har_info[1][1], har_info[1][2])
        if i+1 != len(ene_infos):
            ene_str += ' +'

    return ene_str


# kTP
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


# THERMO
def nasa_polynomial(hform0, hform298, ckin_poly_str):
    """ write the nasa polynomial str
    """
    hf_str = (
        '! Hf(0 K) = {:.2f}, '.format(hform0) +
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


# TRANSPORT
def write_transport_file(ckin_trans_str, ckin_path):
    """ write out the transport file
    """
    if not os.path.exists(ckin_path):
        os.makedirs(ckin_path)
    fpath = os.path.join(ckin_path, 'trans.ckin')
    with open(fpath, 'w') as nasa_file:
        nasa_file.write(ckin_trans_str)
