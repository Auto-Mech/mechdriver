"""
"""


def _write_varecof_inp(varecof_inp, vrc_path):
    """ Write all of the VaReCoF inut files and run the code
    """

    exe_names = ('molpro.sh')

    # Write all of the VaReCoF input files
    for (inp_str, inp_name) in varecof_inp:
        file_name = os.path.join(vrc_path, inp_name)
        # Write the file
        with open(file_name, 'w') as inp_file:
            inp_file.write(inp_str)
        # Make file an executable
        if inp_name in exe_names:
            os.chmod(file_name, mode=os.stat(file_name).st_mode | stat.S_IEXEC)


