""" Handle reading name-label conversion tables for various
    MESS output files
"""

import autoparse.find as apf
import autoparse.pattern as app


# Relabeling functions
def relabel(rxn_ktp_dct, label_dct):
    """ Relabel the rxn ktp dictionaries using the label dictionary
    """

    def _relabel(lbl, label_dct):
        """ Use the label dictionary to change the names in
            MESS output labeling to mech input labeling
        """
        if lbl in label_dct:
            relbl = tuple(label_dct[lbl].split('+'))
        else:
            relbl = (lbl,)
        return relbl

    relab_rxn_ktp_dct = {}
    for rxn, _ktp_dct in rxn_ktp_dct.items():
        rcts, prds, thirdbody = rxn

        relab_rcts = ()
        for rct in rcts:
            relab_rcts += _relabel(rct, label_dct)
        relab_prds = ()
        for prd in prds:
            relab_prds += _relabel(prd, label_dct)

        relab_rxn_ktp_dct[(relab_rcts, relab_prds, thirdbody)] = _ktp_dct

    return relab_rxn_ktp_dct


def name_label_dct(output_str):
    """ Build a dictionary that maps the names of wells and bimols
        barriers provided in the MESS input to the labels that are used
        internally by MESS and printed in the rate output with rate constants.

        Ignores inner and oute barriers, hard to get the block for all of the
        types of output.

        Dictionary is built by parsing the name-label conversion table provided
        in the MESS rate output file and then placing the elements of that
        table into the dictionary.

        :param output_str: string of lines of MESS output file
        :type output_str: str
        :rtype: dict[str: str]
    """

    start_table_ptt = 'Names Translation Tables'
    end_table_ptt = 'Barriers:'
    start_table_ptt2 = 'Well Names Translation'
    end_table_ptt2 = 'End'
    if apf.has_match(start_table_ptt, output_str):
        start_idx, end_idx = None, None
        out_lines = output_str.splitlines()
        for i, line in enumerate(out_lines):
            if start_table_ptt in line:
                start_idx = i
            if end_table_ptt in line:
                end_idx = i
                break

        table_block = '\n'.join(out_lines[start_idx+1: end_idx+1])
        name_ptt = (
            app.capturing(app.one_or_more(app.NONNEWLINE)) +
            app.SPACES +
            '-' +
            app.SPACES +
            app.capturing(app.one_or_more(app.NONNEWLINE))
        )

        name_caps = apf.all_captures(name_ptt, table_block)
        lbl_dct = {lbl.strip(): name.strip() for lbl, name in name_caps}
        
    elif apf.has_match(start_table_ptt2, output_str):
        # older formatting of output table
        start_idx, end_idx = None, None
        out_lines = output_str.splitlines()
        check_outcapture = 0
        for i, line in enumerate(out_lines):
            if start_table_ptt2 in line:
                start_idx = i
            if end_table_ptt2 in line:
                check_outcapture += 1
                if check_outcapture == 2:
                    end_idx = i
                    break

        table_block = '\n'.join(out_lines[start_idx+1: end_idx])
        name_ptt = (
            app.capturing(app.one_or_more(app.NONNEWLINE))
        )
        name_caps = apf.all_captures(name_ptt, table_block)
        lbl_dct = {}
        for entry in name_caps:
            if 'End' not in entry and 'Bimolecular' not in entry:
                lbl, name = entry.strip().split()
                lbl_dct[lbl] = name

    else:
        lbl_dct = None
        print('Warning no name-label table found in output to relabel labels '
              'for rate constants')

    return lbl_dct


def labels(file_str, read_fake=False, mess_file='out'):
    """ Read the labels out of a MESS file
    """

    if mess_file == 'out':
        lbls = _labels_output_string(file_str)
    elif mess_file == 'inp':
        lbls = _labels_input_string(file_str)
    else:
        lbls = ()

    if not read_fake:
        lbls = tuple(lbl for lbl in lbls
                     if 'F' not in lbl and 'f' not in lbl)

    return lbls


def _labels_input_string(inp_str):
    """ Read the labels out of a MESS input file
    """

    def _read_label(line, header):
        """ Get a label from a line
        """
        lbl = None
        idx = 2 if header != 'Barrier' else 4
        line_lst = line.split()
        if len(line_lst) == idx and '!' not in line:
            lbl = line_lst[idx]
        return lbl

    lbls = ()
    for line in inp_str.splitlines():
        if 'Well ' in line:
            lbls += (_read_label(line, 'Well'),)
        elif 'Bimolecular ' in line:
            lbls += (_read_label(line, 'Bimolecular'),)
        elif 'Barrier ' in line:
            lbls += (_read_label(line, 'Barrier'),)

    return lbls


def _labels_output_string(out_str):
    """ Read the labels out of a MESS input file
    """

    lbls = []
    for line in out_str.splitlines():
        if 'T(K)' in line and '->' not in line:
            rxns = line.strip().split()[1:]
            line_lbls = [rxn.split('->') for rxn in rxns]
            line_lbls = [lbl for sublst in line_lbls for lbl in sublst]
            lbls.extend(line_lbls)

    # Remove duplicates and make lst a tuple
    lbls = tuple(set(lbls))

    return lbls
