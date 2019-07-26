""" test autofile.fs
"""
import os
import tempfile
import autofile.fs

PREFIX = tempfile.mkdtemp()
print(PREFIX)


def test__direction():
    """ test autofile.fs.direction
    """
    prefix = os.path.join(PREFIX, 'run_trunk')
    os.mkdir(prefix)

    dir_fs = autofile.fs.direction(prefix)


if __name__ == '__main__':
    test__direction()
