""" fix names in a thermo file with combined names
"""

import os
import sys
import ioformat
import mechanalyzer.builder

PFX = ('s1', 's2', 's3', 's4',
       's5', 's6', 's7', 's8',
       's9', 's10', 's11', 's12')


def main(therm_str, workdir, therm_out_file_name):
    """ Execute
    """

    # Rewrite file string
    new_therm_str = change_names_in_file(therm_str)

    # Write new file
    ioformat.pathtools.write_file(
        new_therm_str, workdir, therm_out_file_name)


def change_names_in_file(therm_str):
    """ amend the names to drop solo- and combined- labels
    """

    def ste_name_line(line):
        if '200.0    1000.0  3000.0      1' in line:
            tmp = line.strip().split()[0]
            tmp2 = tmp.split('-')[1:]
            ste_name = '-'.join(tmp2)
        else:
            ste_name = None
        return ste_name

    def format_name(name):
        name = mechanalyzer.builder.remove_stereo_name_suffix(name)
        name = f'{name:<24s}'
        return name

    new_lines = []
    for line in therm_str.splitlines():
        ste_name = ste_name_line(line)
        if ste_name is not None:
            new_line = format_name(ste_name) + line[24:]
            new_lines.append(new_line)
        else:
            new_lines.append(line)

    return '\n'.join(new_lines)


if __name__ == '__main__':

    # Read file
    CWD = os.getcwd()
    THERM_FILE_NAME = sys.argv[1]
    THERM_OUT_FILE_NAME = THERM_FILE_NAME + '_mod'
    THERM_STR = ioformat.pathtools.read_file(CWD, THERM_FILE_NAME)

    # Run code
    main(THERM_STR, CWD, THERM_OUT_FILE_NAME)
