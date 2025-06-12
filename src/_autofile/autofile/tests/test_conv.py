""" test
"""

import os
import tempfile
from autofile import directory_to_dictionary


PREFIX = tempfile.mkdtemp()
print(PREFIX)


REF_DCT = {
    'A': {
        'type': 'directory',
        'content': {
            'i': {
                'type': 'directory',
                'content': {
                    'q': {
                        'type': 'directory',
                        'content': {
                            'Aiq.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'p': {
                        'type': 'directory',
                        'content': {
                            'Aip.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'Ai.dat': {
                        'type': 'file',
                        'content': '<str>'}}},
            'j': {
                'type': 'directory',
                'content': {
                    'q': {
                        'type': 'directory',
                        'content': {
                            'Ajq.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'Aj.dat': {
                        'type': 'file',
                        'content': '<str>'},
                    'p': {
                        'type': 'directory',
                        'content': {
                            'Ajp.dat': {
                                'type': 'file',
                                'content': '<str>'}}}}}}},
    'B': {
        'type': 'directory',
        'content': {
            'i': {
                'type': 'directory',
                'content': {
                    'q': {
                        'type': 'directory',
                        'content': {
                            'Biq.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'p': {
                        'type': 'directory',
                        'content': {
                            'Bip.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'Bi.dat': {
                        'type': 'file',
                        'content': '<str>'}}},
            'j': {
                'type': 'directory',
                'content': {
                    'q': {
                        'type': 'directory',
                        'content': {
                            'Bjq.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'p': {
                        'type': 'directory',
                        'content': {
                            'Bjp.dat': {
                                'type': 'file',
                                'content': '<str>'}}},
                    'Bj.dat': {'type': 'file', 'content': '<str>'}}}}}
}


def test__():
    """ test directory_to_dictionary
    """

    _build_fs(PREFIX)
    dct = directory_to_dictionary(PREFIX)
    assert REF_DCT == dct


def _build_fs(dir_path):
    """ make a fake filesystem
    """

    locs1 = ('A', 'B')
    locs2 = ('i', 'j')
    locs3 = ('p', 'q')

    for loc1 in locs1:
        path1 = os.path.join(dir_path, loc1)
        os.mkdir(path1)
        for loc2 in locs2:
            path2 = os.path.join(path1, loc2)
            os.mkdir(path2)
            fname2 = os.path.join(
                path2, '{}{}.dat'.format(loc1, loc2))
            with open(fname2, 'w') as fobj:
                fobj.write('<str>')
            for loc3 in locs3:
                path3 = os.path.join(path2, loc3)
                os.mkdir(path3)
                fname3 = os.path.join(
                    path3, '{}{}{}.dat'.format(loc1, loc2, loc3))
                with open(fname3, 'w') as fobj:
                    fobj.write('<str>')


if __name__ == '__main__':
    test__()
