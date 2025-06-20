""" A utility file for various helper functions in this CLI module 
"""

import yaml


def read_mechs_yaml(fname: str, comparison_type: str):
    """ Read the yaml file that specifies all files for each mechanism in a 
        mechanism comparison (either rates or thermo). Also, ensure all files
        are given.

        comparison_type is either 'rates' or 'thermo'
    """ 
    def _mech_is_fully_specified(mech_data, comparison_type):
        """ Ensure all files are given in yaml for mechanism. (Required files
            depends on the comparison type)
        """
        if comparison_type == 'rates':
            ftypes = ['rate_file', 'therm_file', 'species_csv']
        else:
            ftypes = ['therm_file', 'species_csv']
        _fully_specified = True
        if not all(ftype in mech_data for ftype in ftypes):
            _fully_specified = False
            print(
                f'Missing in yaml file {fname}: ',  
                ','.join([ftype for ftype in ftypes if ftype not in mech_data]))
        return _fully_specified

    with open(fname, 'r') as file:
        mechs = yaml.safe_load(file)
    mechs =  {
        key: mech_files for key, mech_files in mechs.items() 
        if _mech_is_fully_specified(mech_files, comparison_type)}

    labels = []
    mech_files = []
    therm_files = []
    csv_files = []
    for key, mech in mechs.items():
        labels.append(key)
        if comparison_type == 'rates':
            mech_files.append(mech['rate_file'])
        therm_files.append(mech['therm_file'])
        csv_files.append(mech['species_csv'])

    return labels, mech_files, therm_files, csv_files
