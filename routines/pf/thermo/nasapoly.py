"""
 Generates NASA Polynomial from MESS+THERMP+PAC99 outputs
"""

def build_polynomial():
    """
    """
   
    print('Generating NASA polynomials at path: {}', nasa_path)

    # Go to NASA path
    pfrunner.go_to_path(nasa_path)

    # Write and run the thermp file to get Hf0k and ...
    pfrunner.write_thermp_inp(spc_dct[spc_name], temps)

    # use new reader function and just pass the 298 k val out of dct
    hf298k = pfrunner.run_thermp(pf_path, nasa_path)
    spc_dct[spc_name]['Hfs'].append(hf298k)

    # Run PAC99 to get a NASA polynomial string in its format
    pac99_str = pfrunner.run_pac(spc_dct[spc_name], nasa_path)
    poly_str = thmroutines.nasapoly.get_pac99_polynomial(pac99_str)

    # Convert the polynomial from PAC99 to CHEMKIN
    chemkin_poly_str = writer.chemkin.run_ckin_poly(
        spc_name, spc_dct[spc_name], poly_str)

    # Write a string for a single spc file to a set
    chemkin_spc_str = chemkin_header_str + chemkin_poly_str
    chemkin_set_str += chemkin_poly_str
            
    # Go back to starting path
    pfrunner.go_to_path(starting_path)


def print_nasa_temps(temps):
    """ Print the polynomial fit temps
    """
    print('Attempting to fit NASA polynomials from',
          '200-1000 and 1000-3000 K ranges using\n',
          'Temps from MESSPF file = {}.'.format(
              ' '.join(('{:.2f}'.format(x) for x in temps))))
