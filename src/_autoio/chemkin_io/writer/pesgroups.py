""" derive pesgroups string for inp writing for mechdriver
"""
# default values

BF_THRESH = 0.1

def write_pes_groups(grp_dct_list, modeltype='rovib_dos', bf_thresh=BF_THRESH):
    """ writes the pes groups file from list of dictionaries            
        :param grp_dct_list: [{'grp': N, 'idxs': [1:1, 1:2, ...], 'peds': [[],[A+B=C+D],...], 'hot': [[C],[]]}, ...]
        :type grps: list(dct)
        :return pes_groups_str: string to write
        :rtype: str
    """
    pg_str = ''
  
    for grp in grp_dct_list:
        pg_str += 'grp {} \n'.format(grp['grp'])
        pg_str += '\t idxs = {} \n'.format(grp['idxs']).replace("'","")
        pg_str += '\t peds = {} \n'.format(grp['peds'])
        pg_str += '\t hot = {} \n'.format(grp['hot'])
        
        if 'modeltype' in grp.keys():
            model = grp['modeltype']
        else:
            model = modeltype
        pg_str += '\t modeltype = {} \n'.format(model).replace("'","")
        pg_str += '\t bf_threshold = {} \n'.format(bf_thresh)
        pg_str += 'end grp \n\n'

    return pg_str