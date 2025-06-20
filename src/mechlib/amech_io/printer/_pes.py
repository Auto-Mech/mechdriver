"""
PES printing
"""

from mechlib.amech_io.printer._print import message


def pes(pes_idx, fml, sub_pes_idx):
    """ a
    """
    message(f'Running PES {pes_idx}: {fml}, SUB PES {sub_pes_idx}',
            newline=1)


def channel(chn_idx, reacs, prods):
    """ a
    """
    rct_str, prd_str = '+'.join(reacs), '+'.join(prods)
    message(f'Running Channel {chn_idx}: {rct_str} = {prd_str}')
