""" Converts the file structure in a given directory to a Python dictionary
"""

from pathlib import Path


def directory_to_dictionary(dir_path):
    """ Build dictionary that maps the branching structure of a Linux
        filesystem into a set of a subdictionaries.

        :param dir_path: root directory path with all files and directories
        :type dir_path: str
        :rtype: dict[str: dict)
    """

    _dic = {}
    for path in Path(dir_path).glob('*'):
        fullpath = path.absolute().resolve()
        if path.is_file():
            with open(fullpath, mode='r', encoding='utf-8') as fobj:
                tmp = fobj.read()
            _dic2 = {
                'type': 'file',
                'content': tmp
            }
        elif path.is_dir():
            _dic2 = {
                'type': 'directory',
                'content': directory_to_dictionary(fullpath)
            }
        _dic[path.name] = _dic2

    return _dic
