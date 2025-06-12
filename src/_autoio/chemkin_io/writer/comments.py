""" Handles the writing of comment strings
"""

from statistics import mean


def get_rxn_cmts_dct(rxn_err_dct=None, rxn_ref_dct=None, rxn_sort_dct=None,
                     rxn_block_cmt=None):
    """ Creates a dict with comments for all rxns based on information about
        the fitting errors, references, and/or sorting. Comments may be
        headers, footers, or inline.

        :param rxn_err_dct: fitting errors for each reaction
        :type rxn_err_dct: dict {rxn: err_dct}
        :param rxn_ref_dct: unknown
        :type rxn_ref_dct:
        :param rxn_sort_dct: unknown
        :type rxn_sort_dct:
        :return rxn_cmts_dct: comments for each reaction
        :rtype: dict {rxn: cmts_dct}
    """

    rxn_err_dct = rxn_err_dct or {}  # set to empty dictionary if None
    rxn_ref_dct = rxn_ref_dct or {}
    rxn_sort_dct = rxn_sort_dct or {}

    # Get all unique reaction keys
    unique_rxns = set(list(rxn_err_dct.keys()) + list(rxn_ref_dct.keys()) +
                      list(rxn_sort_dct.keys()))

    rxn_cmts_dct = {}
    for rxn in unique_rxns:
        cmts_dct = {}

        # Load the data to write the strings
        err_dct = rxn_err_dct.get(rxn)
        sort_dct = rxn_sort_dct.get(rxn)
        # add something here later for references

        # Get the comment strings
        header_cmt = get_header_cmt(sort_dct)
        inline_cmt = get_inline_cmt(sort_dct)
        footer_cmt = get_footer_cmt(err_dct=err_dct)

        # Store the comments
        cmts_dct['header'] = header_cmt
        cmts_dct['inline'] = inline_cmt
        cmts_dct['footer'] = footer_cmt
        rxn_cmts_dct[rxn] = cmts_dct

    # Add the block comment if it was given
    if rxn_block_cmt is not None:
        rxn_cmts_dct['block'] = rxn_block_cmt

    return rxn_cmts_dct


def get_header_cmt(sort_dct=None):
    """ Gets the header comments for a single reaction

        :param sort_dct: sorting comments for a single reaction
        :type sort_dct: dict {'cmts_top': cmt, 'cmts_inline': cmt}
        :return header_cmt: header comment for a single reaction
        :rtype: str
    """

    header_cmt = ''
    if sort_dct is not None:
        header_cmt += sort_dct['cmts_top']

    return header_cmt


def get_inline_cmt(sort_dct=None):
    """ Gets the inline comments for a single reaction

        :param sort_dct: sorting comments for a single reaction
        :type sort_dct: dict {'cmts_top': cmt, 'cmts_inline': cmt}
        :return inline_cmt: inline comment for a single reaction
        :rtype: str
    """

    inline_cmt = ''
    if sort_dct is not None:
        inline_cmt += sort_dct['cmts_inline']

    return inline_cmt


def get_footer_cmt(err_dct=None, ref_dct=None):
    """ Gets the footer comments for a single reaction

        :param err_dct: fitting errors for a single reaction
        :type err_dct: dict {pressure: (temps, errs)}
        :param ref_dct: unknown
        :type ref_dct:
        :return footer_cmt: footer comment for a single reaction
        :rtype: str
    """

    footer_cmt = ''
    if ref_dct is not None:
        pass  # will add this functionality later
    if err_dct is not None:
        footer_cmt += get_err_str(err_dct)

    return footer_cmt


def get_err_str(err_dct):
    """ Gets the string describing fitting errors and ranges for a single reaction

        :param err_dct: fitting error information for a single reaction
        :type err_dct: dict {pressure: (temps, errs)}
        :return err_str: string describing fitting errors for a single reaction
        :rtype: str
    """

    err_str = '! Fitting errors and ranges:\n'

    for pressure, (temps, errs) in err_dct.items():
        max_temp = max(temps)
        min_temp = min(temps)
        mean_err = mean(abs(errs))
        max_err = max(abs(errs))

        if pressure == 'high':
            pressure_str = 'HPL'
        else:
            pressure_str = f'{pressure:.2e} atm'

        err_str += (f'! {pressure_str}: fit betw. {int(min_temp)} and'
                    f' {int(max_temp)} K, MeanAbsErr of {mean_err:.1f}%,'
                    f' MaxAbsErr of {max_err:.1f}%\n')

    return err_str
