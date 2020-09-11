""" Driver for energy transfer evaluations including determination of
    Lennard-Jones parameters
"""

from routines.etrans import lj as ljroutines


def run(spc_dct,
        thy_dct, rxn_lst
        run_inp_dct,
        run_onedmin=True
        write_transport=True):
    """ main driver for etransfer run
    """

    # Pull stuff from dcts for now
    save_prefix = run_inp_dct['save_prefix']
    run_prefix = run_inp_dct['run_prefix']

    # Get the info
    # Read the targets and baths from the input
    TARGET_DCTS = py1dmin.ljparser.read_targets(INPUT_STRING)
    BATH_LST = py1dmin.ljparser.read_baths(INPUT_STRING)

    # Read the run parameters
    THEORY_LEVEL = py1dmin.ljparser.read_theory_level(INPUT_STRING)
    POTENTIAL = py1dmin.ljparser.read_potential(INPUT_STRING)
    NJOBS = py1dmin.ljparser.read_njobs(INPUT_STRING)
    NSAMPS = py1dmin.ljparser.read_nsamps(INPUT_STRING)
    SMIN = py1dmin.ljparser.read_smin(INPUT_STRING)
    SMAX = py1dmin.ljparser.read_smax(INPUT_STRING)
    CONF = py1dmin.ljparser.read_conf(INPUT_STRING)
    TARGET_SAVE_PREFIX = py1dmin.ljparser.read_save_prefix(INPUT_STRING)
    BATH_SAVE_PREFIX = py1dmin.ljparser.read_save_prefix(INPUT_STRING)
    RUN_PREFIX = py1dmin.ljparser.read_run_prefix(INPUT_STRING)

    # Read the theory.dat file into a string
    with open(os.path.join(DRIVE_PATH, 'theory.dat'), 'r') as theoryfile:
        THEORY_STRING = theoryfile.read()

    # Check theory level parameters
    py1dmin.ljparser.check_defined_theory_level_keywords(
        THEORY_STRING, THEORY_LEVEL)

    # Set additional pieces of info using input keywords
    THEORY_INFO = [RUN_PROG, RUN_METHOD, RUN_BASIS, RUN_ORB_REST]

    # Build a list of the species to calculate thermochem for loops below
    spc_queue = parser.species.build_queue(rxn_lst)
    print(spc_queue)

    # Loop over Tasks
    for tsk_lst in trans_tsk_lst:
        # Unpack the options
        [obj, tsk, etrans_keyword_dct] = tsk_lst

        if tsk == 'run_onedmin':
            print(('\n\n------------------------------------------------' +
                   '--------------------------------------'))
            print('\nObtaining LJ-Params using OneDMin')
            for spc_name, _ in spc_queue:
                ljroutines.run_onedmin(
                    spc_name, spc_dct, thy_dct,
                    etrans_keyword_dct)

        if tsk == 'write_transport':
            print(('\n\n------------------------------------------------' +
                   '--------------------------------------'))
            print('\nWriting the CHEMKIN transport file')
            ckin_trans_str = transroutines.write_file()
            writer.ckin.write_transport_file(ckin_trans_str, ckin_path)
