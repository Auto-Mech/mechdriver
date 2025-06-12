""" Parse mechanism file
"""

import os
import mechanalyzer
import ioformat


# Set Paths to test/data directory and output directory
CWD = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(CWD, 'data')

BLD_INP_STR = ioformat.pathtools.read_file(
    DATA_PATH, 'build.dat',
    remove_comments='#',
    remove_whitespace=True)


def test__input_parse():
    """ test mechanalyzer.parser._bld.build_input_string
    """

    ref_inp_dct = {
        'inp_spc': 'species.csv',
        'out_spc': 'species.csv2',
        'inp_mech': 'mechanism.dat',
        'out_mech': 'mechanism.dat2',
        'sort': 'sort.dat',
        'headers': ('smiles', 'inchi', 'mult', 'charge')
    }
    ref_rseries = (
        (('C4H9(1)', 'C4H9(2)'),
         ('O2',),
         None,
         ('addition', 'hydrogen_migration', 'beta_scission')),)

    inp_dct, rseries = mechanalyzer.parser.build_input_file(BLD_INP_STR)

    assert ref_inp_dct == inp_dct
    assert ref_rseries == rseries


def test__reaction_build():
    """ test mechanalyzer.builder.rxn
    """

    spc_dct = {
        'C4H10': {'inchi': 'InChI=1S/C4H10/c1-3-4-2/h3-4H2,1-2H3'},
        'C4H9(1)': {'inchi': 'InChI=1S/C4H9/c1-3-4-2/h1,3-4H2,2H3'},
        'O2': {'inchi': 'InChI=1S/O2/c1-2'},
        'H2': {'inchi': 'InChI=1S/H2/h1H'},
        'H': {'inchi': 'InChI=1S/H'}
    }

    rxn_param_dct = {
        (('C4H10', 'H'), ('C4H9(1)', 'H2'), (None,)): (
            ([1, 0, 0], None, None, None, None, None),)
    }

    rseries = (
        (('C4H9(1)',),
         ('O2',),
         None,
         ('addition', 'hydrogen migration', 'beta scission')),)

    bld_spc_dct, bld_rxn_param_dct = mechanalyzer.builder.rxn.build_mechanism(
        spc_dct, rxn_param_dct, rxn_series=rseries)

    assert set(bld_spc_dct.keys()) == {
        'C4H10',
        'C4H9(1)',
        'O2',
        'H2',
        'H',
        'C4H9(2)',
        'C4H8',
        'C4H9O2',
        'C2H5',
        'C2H4'
    }
    assert set(bld_rxn_param_dct.keys()) == {
        (('C4H10', 'H'), ('C4H9(1)', 'H2'), (None,)),
        (('C4H9(1)',), ('C4H9(1)',), (None,)),
        (('C4H9(1)',), ('C4H9(2)',), (None,)),
        (('C4H9(1)',), ('C4H8', 'H'), (None,)),
        (('C4H9(1)',), ('C2H5', 'C2H4'), (None,)),
        (('C4H9(1)', 'O2'), ('C4H9O2',), (None,))
    }


if __name__ == '__main__':
    test__input_parse()
    test__reaction_build()
