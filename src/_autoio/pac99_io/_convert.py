"""
  Convert the PAC99 polynomial
"""


def pac2ckin_poly(name, atom_dct, pac99_poly_str):
    """ Convert a NASA polynomial from its format given in PAC99 input
        to the format used in ChemKin mechanisms.

        Currently supports C,H,N,P,O,S,F,Cl,Br,I.

        :param name: species name
        :type name: str
        :param atom_dct: number of each atom in the species
        :type atom_dct: dict[str: int]
        :param pac99_poly_str: PAC99-format NASA polynomial
        :type pac99_poly_str: str
        :rtype: str
    """

    # Build the elemental string and determine if polynomial can be converted
    elem_str = _write_composition_str(atom_dct)

    # If possible, parse the remaining output and write the converted string
    if elem_str is not None:
        las, has = _parse_coefficients(pac99_poly_str)
        lowt, breakt, hight = _parse_temperatures(pac99_poly_str)

        line1 = (
            f'{name.strip():<24s}' +      # name string (assume 9-char buffer)
            # '  ' +                      # 2-char buffer
            # '      ' +                  # unused 6-char comment
            f'{elem_str:<20s}' +
            'G'
            f'{lowt:>10.1f}{hight:>10.1f}{breakt:>8.1f}'
            '      ' +                  # unused 6-char comment
            '1\n'
        )
        line2 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    2\n" % (
            has[0], has[1], has[2], has[3], has[4])
        line3 = "% 15.8E% 15.8E% 15.8E% 15.8E% 15.8E    3\n" % (
            has[5], has[6], las[0], las[1], las[2])
        line4 = "% 15.8E% 15.8E% 15.8E% 15.8E                   4\n" % (
            las[3], las[4], las[5], las[6])
        full_line = line1 + line2 + line3 + line4
    else:
        full_line = None

    return full_line


# Specific writers
def _write_composition_str(atom_dct):
    """ Build the string that details the elemental composition of the species
        Currently supports C,H,N,P,O,S,F,Cl,Br,I.

        :param atom_dct: number of each atom in the species
        :type atom_dct: dict[str: int]
        :rtype: str
    """

    # Build the string that details the elemental composition of the species
    num_c = atom_dct.get('C', 0)
    num_h = atom_dct.get('H', 0)
    num_n = atom_dct.get('N', 0)
    num_p = atom_dct.get('P', 0)
    num_o = atom_dct.get('O', 0)
    num_s = atom_dct.get('S', 0)
    num_f = atom_dct.get('F', 0)
    num_cl = atom_dct.get('Cl', 0)
    num_br = atom_dct.get('Br', 0)
    num_i = atom_dct.get('I', 0)

    elem_str = ''
    if num_c > 0:
        elem_str += f'C{num_c:>4d}'
    if num_h > 0:
        elem_str += f'H{num_h:>4d}'
    if num_n > 0:
        elem_str += f'N{num_n:>4d}'
    if num_p > 0:
        elem_str += f'P{num_p:>4d}'
    if num_o > 0:
        elem_str += f'O{num_o:>4d}'
    if num_s > 0:
        elem_str += f'S{num_s:>4d}'
    if num_f > 0:
        elem_str += f'F{num_f:>4d}'
    if num_cl > 0:
        elem_str += f'CL{num_cl:>3d}'
    if num_br > 0:
        elem_str += f'BR{num_br:>3d}'
    if num_i > 0:
        elem_str += f'I{num_i:>4d}'

    elem_str = elem_str.strip()

    if len(elem_str) > 20:
        print('WARNING: TOO MANY ELEMS TO WRITE COMPOSITION STRING OF POLYNOM')
        elem_str = None

    return elem_str


# Parsers
def _parse_temperatures(pac_poly_str):
    """ Parse the temperatures that define the low- and high-temperature range
        of the polynomial in the .c97 output file.

        Returns the low-temp, break-temp, and high-temp

         :param pac_poly_str: string polynomial
         :type string: str
         :rtype: tuple(tuple(float))
    """

    # Parse the lines of the pac string containing the desired coefficients
    pac_lines = pac_poly_str.splitlines()
    line1 = pac_lines[5].strip().split()  # contains low and break temp
    line2 = pac_lines[8].strip().split()  # contains break and high temp

    return float(line1[0]), float(line1[1]), float(line2[1])


def _parse_coefficients(pac_poly_str):
    """  Parse the coefficients from teh .c97 output file

         :param out_str: string polynomial
         :type string: str
         :rtype: tuple(tuple(float))
    """

    def _parse_line16(string):
        """ Parse the values of a string containing exponental numbers
            of length 16 chars, where numbers may not have a space in between.
        """

        assert len(string) % 16 == 0, 'String should have 16n chararacters'

        # Build a list of values from the string
        string2 = string.replace('D', 'E')
        nchunks = len(string2) // 16
        vals = [0.0 for i in range(nchunks)]
        for i in range(nchunks):
            vals[i] = float(string2[i*16: (i+1)*16])

        return vals

    # Parse the lines of the pac string containing the desired coefficients
    pac_lines = pac_poly_str.splitlines()

    las = [0.0 for i in range(7)]
    has = [0.0 for i in range(7)]
    las[0:5] = _parse_line16(pac_lines[6][0:80])
    las[5:7] = _parse_line16(pac_lines[7][48:80])
    has[0:5] = _parse_line16(pac_lines[9][0:80])
    has[5:7] = _parse_line16(pac_lines[10][48:80])

    return tuple(las), tuple(has)
