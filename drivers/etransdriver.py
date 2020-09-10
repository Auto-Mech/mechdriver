""" Driver for energy transfer evaluations including determination of
    Lennard-Jones parameters
"""

# from routines.etrans import x
# from routines.etrans import etransrunner


def run(spc_dct,
        thy_dct, rxn_lst
        run_inp_dct):
    """ main driver for etransfer run
    """

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']
    run_prefix = run_inp_dct['run_prefix']

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = parser.species.build_queue(rxn_lst)
    print(spc_queue)
   
    # Loop over Tasks
    print('\nRunning energy transfer tasks given in the input...')
    for tsk_lst in es_tsk_lst:

        # Unpack the options
        [obj, tsk, etrans_keyword_dct] = tsk_lst

        spc_queue = parser.species.build_spc_queue(rxn_lst)
