"""
  Read the output of the PAC99
"""


def nasa_polynomial(output_str):
    """ Parse the NASA polynomial from the PAC99 output file.

        :param output_str: string for the output file
        :type output_str: str
        :rtype: str
    """

    lines = output_str.splitlines()

    return '\n'.join([lines[i] for i in range(11)])
