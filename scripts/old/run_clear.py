""" script
"""
import os
import autofile


def clear_running_status(path):
    """ clear running status
    """
    inf_paths = find_files('run.yaml', path=path)
    for inf_path in inf_paths:
        inf_str = autofile.file.read_file(inf_path)
        inf_obj = autofile.file.read.information(inf_str)
        if inf_obj.status == 'running':
            inf_obj.status = 'failed'
        inf_str = autofile.file.write.information(inf_obj)
        autofile.file.write_file(inf_path, inf_str)


def find_files(file_name, path):
    """ file files by name
    """
    file_paths = []
    for dir_path, _, file_names in os.walk(path):
        if file_name in file_names:
            file_paths.append(os.path.join(dir_path, file_name))
    return tuple(file_paths)


if __name__ == '__main__':
    clear_running_status('.')

