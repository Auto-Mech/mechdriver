"""
Executes the automation part of 1DMin
"""

import os
import sys
import moldr
import autofile
import py1dmin.ljroutines
import py1dmin.ljparser


# Set paths to various working directories
DRIVE_PATH = os.getcwd()
GEOM_PATH = os.path.join(DRIVE_PATH, 'geoms')

# Read the input file into a string
with open(os.path.join(DRIVE_PATH, 'input.dat'), 'r') as infile:
    INPUT_STRING = infile.read()

# Check to see if the parameters are defined
py1dmin.ljparser.check_defined_sections(INPUT_STRING)
py1dmin.ljparser.check_defined_lennard_jones_keywords(INPUT_STRING)

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

# Read the theory level parameters
RUN_PROG = py1dmin.ljparser.read_program(THEORY_STRING, THEORY_LEVEL)
RUN_METHOD = py1dmin.ljparser.read_method(THEORY_STRING, THEORY_LEVEL)
RUN_BASIS = py1dmin.ljparser.read_basis(THEORY_STRING, THEORY_LEVEL)
RUN_ORB_REST = py1dmin.ljparser.read_orb_restrict(THEORY_STRING, THEORY_LEVEL)
RUN_MEMORY = py1dmin.ljparser.read_memory(THEORY_STRING, THEORY_LEVEL)
RUN_NPROCS = py1dmin.ljparser.read_nprocs(THEORY_STRING, THEORY_LEVEL)

# Set additional pieces of info using input keywords
THEORY_INFO = [RUN_PROG, RUN_METHOD, RUN_BASIS, RUN_ORB_REST]


# Loop over the species and launch all of the desired jobs
for target, target_info in TARGET_DCTS.items():

    # Set the new information
    RUN_CHG = target_info[1] + BATH_LST[1]
    RUN_MLT = max([target_info[2], BATH_LST[2]])
    BATH_INFO = BATH_LST

    # Search save file system for LJ params
    SIGMA, EPSILON = py1dmin.ljroutines.read_lj_from_save(
        TARGET_SAVE_PREFIX, target_info, THEORY_INFO)

    # Run 1DMin to calculate the LJ parameters
    if SIGMA is None and EPSILON is None:  # or NSAMP <= nsamp_run:

        # Obtain the geometry for the target and bath
        TARGET_GEO = py1dmin.ljroutines.get_geometry(
            target_info,
            THEORY_INFO,
            TARGET_SAVE_PREFIX,
            geom_dct=moldr.util.geometry_dictionary(GEOM_PATH),
            conf=CONF)
        BATH_GEO = py1dmin.ljroutines.get_geometry(
            BATH_INFO,
            THEORY_INFO,
            BATH_SAVE_PREFIX,
            geom_dct=moldr.util.geometry_dictionary(GEOM_PATH),
            conf=CONF)

        # Write the params to the run file system
        FS_THEORY_INFO = [THEORY_INFO[1],
                          THEORY_INFO[2],
                          moldr.util.orbital_restriction(
                              target_info, THEORY_INFO)]
        tgt_run_fs = autofile.fs.species(RUN_PREFIX)
        tgt_run_fs.leaf.create(target_info)
        tgt_run_path = tgt_run_fs.leaf.path(target_info)
        etrans_run_fs = autofile.fs.energy_transfer(tgt_run_path)
        etrans_run_path = etrans_run_fs.leaf.path(FS_THEORY_INFO)
        etrans_run_fs.leaf.create(FS_THEORY_INFO)

        # Run an instancw of 1DMin for each processor
        for i in range(NJOBS):

            # Build run directory
            job_dir_path = os.path.join(
                etrans_run_path, 'run{0}'.format(str(i+1)))
            os.mkdir(job_dir_path)
            print('\n\nWriting files to'+job_dir_path)

            # Write the 1DMin input file
            print('  Writing input files...')
            py1dmin.ljroutines.write_input(
                job_dir_path, NSAMPS,
                target_name='target.xyz', bath_name='bath.xyz',
                smin=SMIN, smax=SMAX)

            # Write the geometry files
            print('  Writing xyz files for target and bath....')
            py1dmin.ljroutines.write_xyz(
                job_dir_path, TARGET_GEO, BATH_GEO)

            # Write the electronic structure template file
            print('  Writing electronic structure submission inp template...')
            py1dmin.ljroutines.write_elstruct_inp(
                job_dir_path,
                RUN_CHG, RUN_MLT, RUN_METHOD, RUN_BASIS, THEORY_INFO,
                RUN_PROG, RUN_MEMORY)

            # Write the electronic structure sumbission script
            print('  Writing electronic structure submission script...')
            py1dmin.ljroutines.write_elstruct_sub(
                job_dir_path, DRIVE_PATH, RUN_PROG)

        # Submit the job
        print('\n\nRunning each OneDMin job...')
        py1dmin.ljroutines.submit_job(DRIVE_PATH, etrans_run_path, NJOBS)

        # Read the lj from all the outputs
        print('\n\nAll OneDMin jobs finished.')
        print('\nReading the Lennard-Jones parameters...')
        SIGMA, EPSILONS = py1dmin.ljroutines.lj_parameters(etrans_run_path)
        if SIGMA is None and EPSILON is None:
            print('\nNo Lennard-Jones parameters found.')
            print('\n\nExiting OneDMin...')
            sys.exit()
        else:
            # Grab the geometries and zero-energies
            print('\nReading the Lennard-Jones Potential Well Geometries...')
            GEOMS = py1dmin.ljroutines.lj_well_geometries(
                etrans_run_path)
            print('\nReading the Zero-Energies...')
            ENES = py1dmin.ljroutines.zero_energies(
                etrans_run_path)

        # Write the params to the save file system
        py1dmin.ljroutines.write_lj_to_save(
            SIGMA, EPSILON,
            TARGET_SAVE_PREFIX, target_info, THEORY_INFO)
        print('\n\nExiting OneDMin...')
