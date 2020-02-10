"""
Analyze the input mechanism and output a new one to be used by the program
"""

def write_stereo_csv_file(spc_dct):
    """ After checking stereochemistry, write a new CSV file with stero inchis
    """
    # Get new stereochemically accurate species

    # Write the new CSV string 
    csv_str = ''
    return csv_str


def write_mech_file(pes_dct):
    """ Write a new CHEMKIN Mechanims file with a PES dictionary
    """
    mech_str = ''
    # call a similar function to the print_pes_channels
    for idx, formula in enumerate(pes_dct):
        # Print the header for each PES (by formula)
        mech_str += '! PES {}: {}\n'.format(idx+1, formula)
        mech_str += reaction_string(pes_dct[formula])    

    return mech_str


def reaction_string(pes):
    """ Get a string from reactants and products (maybe move to chemkin_io)
    """
    pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']
    pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
    pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
    chn_str = ''
    for chn_idx, _ in enumerate(pes_rxn_name_lst):
        chn_str += 'Channel {}:\n{} = {}     1.000  0.00  0.00'.format(
            chn_idx+1,
            ' + '.join(pes_rct_names_lst[chn_idx]),
            ' + '.join(pes_prd_names_lst[chn_idx])))

    return chn_str




