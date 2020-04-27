"""
  Library for writing output in CHEMKIN formatted files
"""

import os
import automol
import routines.pf.thermo


def get_ckin_ene_lvl_str(pf_levels, ene_coeff):
    """ Write the comment lines for the enrgy lvls for ckin
    """
    ene_thy_info = pf_levels[1]
    ene_ref_thy_info = pf_levels[0]
    ene_strl = []
    ene_str = '! energy level:'
    ene_strl.append(' {:.2f} x {}{}/{}//{}{}/{}\n'.format(
        ene_coeff[0],
        ene_thy_info[3],
        ene_thy_info[1],
        ene_thy_info[2],
        ene_ref_thy_info[3],
        ene_ref_thy_info[1],
        ene_ref_thy_info[2]))
    ene_str += '!               '.join(ene_strl)

    return ene_str


def run_ckin_header(pf_info, spc_model):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """
    tors_model, vib_model, _ = spc_model
    [geo_info, ene_info, har_info,
     vpt2_info, _, tors_info] = pf_info
    sp_str = False
    # [geo_info, ene_info, har_info,
    #  vpt2_info, _, tors_info, sp_str] = pf_info

    # Convert the pac99 polynomial to chemkin polynomial
    chemkin_header_str = '! vib model: {0}\n'.format(vib_model)
    chemkin_header_str += '! tors model: {0}\n'.format(tors_model)
    if har_info:
        chemkin_header_str += '! ene//geo level: {}{}/{}//{}{}/{}\n'.format(
            ene_info[3], ene_info[1], ene_info[2],
            geo_info[3], geo_info[1], geo_info[2])
    if tors_info:
        chemkin_header_str += '! tors level: {}{}/{}//{}{}/{}\n'.format(
            tors_info[1][3], tors_info[1][1], tors_info[1][2],
            tors_info[0][3], tors_info[0][1], tors_info[0][2])
    if vpt2_info:
        chemkin_header_str += '! vpt2 level: {}{}/{}\n'.format(
            vpt2_info[3], vpt2_info[1], vpt2_info[2])
    if sp_str:
        chemkin_header_str += sp_str

    return chemkin_header_str


def run_ckin_poly(spc, spc_dct_i, pac99_poly_str):
    """ prepare chemkin header info and convert pac 99 format to chemkin format
    """

    hf_str = '! Hf(0 K) = {:.2f}, Hf(298 K) = {:.2f} kcal/mol\n'.format(
        float(spc_dct_i['Hfs'][0]), float(spc_dct_i['Hfs'][1]))
    ich = spc_dct_i['ich']
    formula_dct = automol.inchi.formula_dct(ich)
    chemkin_poly_str = routines.pf.thermo.nasapoly.convert_pac_to_chemkin(
        spc, formula_dct, hf_str, pac99_poly_str)
    print('\nCHEMKIN Polynomial:')
    print(chemkin_poly_str)
    return chemkin_poly_str


def write_nasa_file(spc_dct_i, ckin_path, nasa_path, chemkin_poly_str):
    """ write out the nasa polynomials
    """
    ich = spc_dct_i['ich']
    formula = automol.inchi.formula(ich)
    with open(os.path.join(nasa_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
    with open(os.path.join(ckin_path, formula+'.ckin'), 'w') as nasa_file:
        nasa_file.write(chemkin_poly_str)
