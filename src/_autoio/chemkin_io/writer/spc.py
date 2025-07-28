""" Writes the parts of a Chemkin file pertaining to species, i.e., the
    species block and the elements block
"""

from automol.chi import formula
def write_species(mech_spc_dct):
    """ Writes the contents of a mech_spc_dct in Chemkin format, with notes

        :param mech_spc_dct: species data for a mechanism
        :type mech_spc_dct: dct {spc_name:data}
        :return spc_str: str containing the species information
        :rtype: str
    """

    # Get the max species name length
    max_spc_len = 0
    for spc in mech_spc_dct.keys():
        max_spc_len = max(len(spc), max_spc_len)

    # Get the max SMILES name length
    max_smiles_len = 0
    for spc, ident_dct in mech_spc_dct.items():
        max_smiles_len = max(len(ident_dct['smiles']), max_smiles_len)

    buffer = 5

    # Write the spc_str
    spc_str = ''
    for spc, ident_dct in mech_spc_dct.items():
        # Determine what ChI prefix to write ('AMChI' or 'InChI')
        if 'AMChI' in ident_dct['inchi']:
            chi_prefix = 'AMChI: '
        else:
            chi_prefix = 'InChI: '
        spc_str += (
            '{0:<'+str(max_spc_len+buffer)+'s}{1:>9s}{2:<' +
            str(max_smiles_len+buffer)+'s}{3:>11s}{4:<9s}\n').format(
                spc, '! SMILES: ',
                ident_dct['smiles'],
                chi_prefix,
                ident_dct['inchi'])

    return spc_str


def write_elements(mech_spc_dct):
    """ Writes the unique elements in a mech_spc_dct to a string

        :param mech_spc_dct: species data for a mechanism
        :type mech_spc_dct: dct {spc_name:data}
        :return elem_str: a str containing the list of unique elements
        :rtype: str
    """
    def _unique_elems(mech_spc_dct):
        """ Gets the unique elements in a mech_spc_dct
        """

        # Build a non-unique list of all elements in the set
        all_elems = []
        for _, spc_dct in mech_spc_dct.items():
            if 'fml' in spc_dct:
                fml = spc_dct['fml']  # assumes that the fml exists!
            else:
                fml = formula(spc_dct['inchi'])
            all_elems.extend(list(fml.keys()))

        # Get the unique elements via a set; sorted otherwise random order
        unique_elems = sorted(list(set(all_elems)))

        return unique_elems

    unique_elems = _unique_elems(mech_spc_dct)
    elem_str = ''
    for unique_elem in unique_elems:
        elem_str += unique_elem + '\n'

    return elem_str
