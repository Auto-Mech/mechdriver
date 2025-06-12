""" test ioformat.pathtools
"""

import os
import tempfile
import numpy
import ioformat


TMP_DIR = tempfile.mkdtemp()


# Generic/Text file
FILE_NAME = 'tmp_file.dat'
FILE_STR = (
    '  ## Gradient (Symmetry 0) ##\n'
    '  Irrep: 1 Size: 3 x 3\n'
    '\n'
    '            1         2         3\n'
    '\n'
    '    1     0.000     0.000     0.997\n'
    '    2     0.000    -0.749    -0.488\n'
    '    3    -0.000     0.749    -0.488\n\n\n')
FILE2_STR = (
    'Irrep: 1 Size: 3 x 3\n'
    '1         2         3\n'
    '1     0.000     0.000     0.997\n'
    '2     0.000    -0.749    -0.488\n'
    '3    -0.000     0.749    -0.488\n')

# NumPy file
NP_NAME = 'tmp_numpy.dat'
NP_ARR = numpy.array(
    ((-4.0048955763, -0.3439866053, -0.0021431734),
     (-1.3627056155, -0.3412713280, 0.0239463418),
     (-4.7435343957, 1.4733340928, 0.7491098889),
     (-4.7435373042, -1.9674678465, 1.1075144307),
     (-4.6638955748, -0.5501793084, -1.9816675556),
     (-0.8648060003, -0.1539639444, 1.8221471090)))

# JSON file
JSON_NAME = 'tmp_json.dat'
JSON_DCT = {
    'a': 1,
    'b': {
        'ba': 10
    },
    'c': {
        'ca': {
            'caa': 100
        }
    }
}

NONEXIST_FILE_NAME = 'nofile.dat'


def test__go_to():
    """ test ioformat.pathtools.current_path
        test ioformat.pathtools.prepare_path
        test ioformat.pathtools.go_to
    """

    cur_path = ioformat.pathtools.current_path()

    tmp_path = ioformat.pathtools.prepare_path((cur_path, TMP_DIR), make=False)
    assert not os.path.exists(tmp_path)

    tmp_path = ioformat.pathtools.prepare_path((cur_path, TMP_DIR), make=True)
    ioformat.pathtools.go_to(tmp_path)
    assert os.getcwd() == tmp_path
    ioformat.pathtools.go_to(cur_path)
    assert os.getcwd() == cur_path


def test__pathtools():
    """ test ioformat.pathtools.write_file
        test ioformat.pathtools.read_file
    """

    # Read file and return same
    ioformat.pathtools.write_file(FILE_STR, TMP_DIR, FILE_NAME)
    file_str = ioformat.pathtools.read_file(TMP_DIR, FILE_NAME)
    assert file_str == FILE_STR

    # Reread the file with certain components removed
    file2_str = ioformat.pathtools.read_file(
        TMP_DIR, FILE_NAME, remove_comments='#', remove_whitespace=True)
    assert file2_str == FILE2_STR

    # Read a file from a path that does not exist
    file3_str = ioformat.pathtools.read_file(TMP_DIR, NONEXIST_FILE_NAME)
    assert file3_str is None


def test__numpy_file():
    """ test ioformat.pathtools.write_numpy_file
        test ioformat.pathtools.read_numpy_file
    """

    ioformat.pathtools.write_numpy_file(NP_ARR, TMP_DIR, NP_NAME)
    arr = ioformat.pathtools.read_numpy_file(TMP_DIR, NP_NAME)
    assert numpy.allclose(NP_ARR, arr)

    # Read a file from a path that does not exist
    arr2 = ioformat.pathtools.read_numpy_file(TMP_DIR, NONEXIST_FILE_NAME)
    assert arr2 is None


def test__json_file():
    """ test ioformat.pathtools.write_json_file
        test ioformat.pathtools.read_json_file
    """

    ioformat.pathtools.write_json_file(JSON_DCT, TMP_DIR, JSON_NAME)
    dct = ioformat.pathtools.read_json_file(TMP_DIR, JSON_NAME)
    assert JSON_DCT == dct

    # Read a file from a path that does not exist
    dct2 = ioformat.pathtools.read_json_file(TMP_DIR, NONEXIST_FILE_NAME)
    assert dct2 is None
