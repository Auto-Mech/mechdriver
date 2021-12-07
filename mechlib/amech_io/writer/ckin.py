"""
  Library for writing output in CHEMKIN formatted files
"""

import os
from phydat import phycon


# COMMON HEADER STUFF FOR kTP and THERMO CKIN FILES
def model_header(spc_mods, spc_mod_dct, refscheme=''):
    """ Write a model header for multiple models
    """
    mod_str = ''
    for spc_mod in spc_mods:
        mod_str += _model_header(spc_mod_dct[spc_mod], refscheme=refscheme)
    mod_str += '\n\n'

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
    chemkin_header_str = f'! vib model: {vib_model}\n'
    chemkin_header_str += f'! tors model: {tors_model}\n'
    chemkin_header_str += f'! sym model: {sym_model}\n'

    # Energy level string
    if har_info is not None and ene_infos:
        chemkin_header_str += '! energy level:\n'
        for i, ene_inf in enumerate(ene_infos):
            if i == 0:
                chemkin_header_str += '!    ' + f'{ene_inf[1][0]:>6.3f} x '
            else:
                chemkin_header_str += '! +  ' + f'{ene_inf[1][0]:>6.3f} x '
            chemkin_header_str += (
                f'{ene_inf[1][1][3]}{ene_inf[1][1][1]}/{ene_inf[1][1][2]}'
                '//'
                f'{har_info[1][3]}{har_info[1][1]}/{har_info[1][2]}\n'
            )

    # Tors level string
    if tors_model != 'rigid':
        chemkin_header_str += (
            '! tors level: '
            f'{tors_ene_info[1][3]}{tors_ene_info[1][1]}/{tors_ene_info[1][2]}'
            '//'
            f'{tors_geo_info[1][3]}{tors_geo_info[1][1]}/{tors_geo_info[1][2]}'
            '\n'
        )
    # if vpt2_info is not None:
    #     chemkin_header_str += '! vpt2 level: {}/{}\n'.format(
    #         vpt2_info[1][1], vpt2_info[1][2])
    if refscheme:
        chemkin_header_str += f'! CBH reference scheme: {refscheme}\n'

    return chemkin_header_str


# kTP


# THERMO
# combine with nasapoly str
def nasa_polynomial(hform0, hform298, ckin_poly_str):
    """ write the nasa polynomial str
    """
    hf_str = (
        f'! Hf(0 K) = {(hform0 * phycon.EH2KCAL):.2f}, '
        f'! Hf(298 K) = {(hform298 * phycon.EH2KCAL):.2f} kcal/mol\n'
    )
    return hf_str + ckin_poly_str


# prob can handle with autorun or ioformat function func
def write_nasa_file(ckin_nasa_str, ckin_path, idx=None):
    """ write out the nasa polynomials
    """
    if not os.path.exists(ckin_path):
        os.makedirs(ckin_path)
    fpath = os.path.join(ckin_path, 'all_therm.ckin')
    if idx is not None:
        fpath += f'_{idx:g}'
    # Add the REACTION section header
    ckin_str = 'THERMO' + '\n\n\n'
    ckin_str += ckin_nasa_str + '\n\n'
    ckin_str += 'END' + '\n\n'

    with open(fpath, mode='w', encoding='utf-8') as nasa_file:
        nasa_file.write(ckin_str)


# TRANSPORT
# prob can handle with autorun func
def write_transport_file(ckin_trans_str, ckin_path):
    """ write out the transport file
    """
    if not os.path.exists(ckin_path):
        os.makedirs(ckin_path)
    fpath = os.path.join(ckin_path, 'trans.ckin')
    with open(fpath, mode='w', encoding='utf-8') as nasa_file:
        nasa_file.write(ckin_trans_str)
