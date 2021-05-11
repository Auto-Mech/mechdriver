"""
PES printing
"""

from mechlib.amech_io.printer import message



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
