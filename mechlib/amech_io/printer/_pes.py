"""
PES printing
"""

from mechlib.amech_io.printer import message


def pes_channels(pes_dct):
    """ Print the PES
    """

    message('Sorted Mechanism read from file:', newline=1, indent=1)
    for pes_idx, formula in enumerate(pes_dct):
        message('! PES:', pes_idx+1, formula)
        pes_rxn_name_lst = pes_dct[formula]['rxn_name_lst']
        pes_rct_names_lst = pes_dct[formula]['rct_names_lst']
        pes_prd_names_lst = pes_dct[formula]['prd_names_lst']
        for chn_idx, _ in enumerate(pes_rxn_name_lst):
            message('  {} = {}   1.0 0.0 0.0'.format(
                ' + '.join(pes_rct_names_lst[chn_idx]),
                ' + '.join(pes_prd_names_lst[chn_idx])))


def pes(pes_idx, formula, sub_pes_idx):
    """ a
    """
    message('Running PES {}: {}, SUB PES {}'.format(
                  pes_idx, formula, sub_pes_idx), newline=1)

def channel(chn_idx, reacs, prods):
    """ a
    """
    message('Running Channel {}: {} = {}'.format(
        chn_idx, '+'.join(reacs), '+'.join(prods)), indent=1)
