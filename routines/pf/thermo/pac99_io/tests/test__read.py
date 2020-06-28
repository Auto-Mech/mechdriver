""" Test pac99_io.reader
"""

import pac99_io


def _read_file(file_name):
    with open(file_name, encoding='utf8', errors='ignore') as file_obj:
        file_str = file_obj.read()
    return file_str


# Set paths
PATH = os.path.dirname(os.path.realpath(__file__))
DATA_PATH = os.path.join(PATH, 'data')
PAC99_OUT_NAME = 'CH4.o97'

# Read mechanism files
PAC99_OUT_STR = _read_file(
    os.path.join(DATA_PATH, PAC99_OUT_NAME))


def test__polynomial():
    """ test pac99_io.reader.__
    """

    nasa_poly = pac99_io.reader.nasa_polynomial(output_str)

    ref_nasa_poly = """
    """

    assert nasa_poly == ref_nasa_poly


if __name__ == '__main__':
    test__polynomial()
