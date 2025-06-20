""" read and write to files
"""
import os


def read_file(file_path):
    """ read a file as a string

    :param file_path: path of file to be read
    :type file_path: str
    :return: file contents
    :rtype: str
    """
    assert os.path.isfile(file_path)
    with open(file_path, mode='r', encoding='utf-8') as file_obj:
        file_str = file_obj.read()
    return file_str


def write_file(file_path, string):
    """ write a string to a file

    :param file_path: path of file to be written
    :type file_path: str
    :param file_path: string to be written
    :type file_path: str
    """
    with open(file_path, mode='w', encoding='utf-8') as file_obj:
        file_obj.write(string)
