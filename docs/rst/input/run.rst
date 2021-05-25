
run.dat
=======

The central input where the user specifies what will be run by the program
    (1) objects, e.g. species, reactions, PESs
    (2) drivers and tasks
    (3) levels of theoretical treatment.

For many of the above, the user will supply labels that correspond to blocks of info that is pulled from the other input file.

The run.dat file is comprised of four sections:

    (1) input
    (2) pes/spc
    (3) els
    (4) trans
    (5) thermo
    (6) ktp
    (7) process

where each section is specified in run.dat as

Input Block
-----------

Simple block that specifies

Example::

    input
        mech = chemkin
        run_prefix = <path/to/run/prefix>
        save_prefix = <path/to/save/prefix>
    end input


.. csv-table:: keywords for theory blocks
    :file: tables/runinp_keys.csv
    :header-rows: 1
    :widths: 10, 10, 10, 10


Chemistry Blocks
----------------

Tells the code what to do calculations for.

Currently both can be given in an output; however the PES block is given by
priority and the SPC block is combined with it. Certain cases where a block
is required.

Example (put one with both PES and SPC)


Potential Energy Surface Block
------------------------------

Lists the Potential Energy Surfaces and constituent channels
(composed of species and transition states) for which to perform calculations.

The indices for the PESs and channels correspond to the those provided in input
`mechanism.dat` file. Note that multiple PES-Channel combinations can be
given in one input.

Required for:
    kTPDriver

General Format::

    pes
        <PES indices>: <Channel indices>
    end

Example 1::

    pes
        1: 1         # Runs Channel 1 on PES 1
        1-2: 1-2     # Runs Channels 1,2 on PESs 1,2
        1,3: 1-4     # Runs Channels 1,2,34 on PESs 1,3
        1,3: 1-4,6   # Runs Channels 1,2,34 on PESs 1,3
    end

Example 2::

    pes
       1: 1-5
       2: 1,2-10
       5: 8
    end

Species Block
-------------

Similar to the PES Block, this section lists the species for which to perform calculations.

The indices for the species correspond to the order of the species provided in the
input `species.csv` file. Note that species can be extended across multiple lines

Required for:
    ThermoDriver

General Format::

    spc
        <Species indices>
    end

Example::

    spc
        1          # Runs species 1
        1,3        # Runs species 1,3
        1-4        # Runs species 1,2,3,4
        1,3-5      # Runs species 1,3,4,5
        1-3,6      # Runs species 1,2,3,6
        1-3,5-7    # Runs species 1,2,3,5,6,7
    end


Task Block
----------

The overall section is

General format::

    tsks
        <object>  <task>  <keyword1=value  keyword2=value …>
        <object>  <task>  <keyword1=value  keyword2=value …>
        <object>  <task>  <keyword1=value  keyword2=value …>
        …
    end

Above, the <object> is either spc (species, i.e. reacs, prods) or ts

Above the <task> is what electronic structure calculation to be run on object.

keyword=value cannot have spaces in between them.

Each task is given in the following format <obj>_<job>
 
Electronic Structure Driver Task Block
---------------------------------------

Specifies electronic structure tasks

The order of the tasks in the input block corresponds to the order the tasks
are conducted in the order they are calculated.

Any number of tasks can be given in the file.


Explanation of keywords

init: geom
find: ts
conf: pucker, samp, energy, grad, hess, vpt2, prop, opt
hr: scan, grad, hess, energy, vpt2, reopt
tau: samp, energy, grad, hess
rpath: scan, energy, grad, hess:

init: task where no info possibly exists using InChI or input geometry
find: TS searching algorithm
conf: conformer of species/TS defined torsional angles of internal rotors and ring torsions
hr: all geometries along the motion of an internal rotor

energy: calculates electronic energy
grad: calculates molecular gradient in Cartesian coordinates
hess: calculates molecular hessian in Cartesian coordinates
vpt2: calculates VPT2 equations and saves several components
prop: calculates molecular properties including dipole moment and polarizability
opt: optimizes some geometry
reopt: optimizes the minimum on the HR potential if a new min is found
samp: optimizes conformers according Monte Carlo sampling routine
scan: scans across some specified internal coordinate

| key descriptions:
| inplvl = electronic structure method where geometry is read
|          either for a guess (searching for scan) or for a calc on top (energy calc)
| runlvl = electronic structure method to calculate desired data for
| tors_model = model used to construct constraints for scan

for electronic structure method, these should correspond to those that defines blocks in 
theory.dat


kTPDriver and ThermDriver Task Blocks
-------------------------------------

The order given in the task block is not used anywhere. We use our own
internal order: write -> run -> fit

| write_mess: Reads and processes all information according to models specified
|             and then builds the required MESS file
|             give path to where the MESS file is saved in run-filesys
| run_mess: Run the MESS file existing in the filesystem
| run_fits: Fit the rate constants (thermochem) to functional forms and output
|           mechanism files with this data

kin_model: conditions model (kin mod) to use that is defined in models.dat
spc_model: conditions model (spc mod) to use that is defined in models.dat
overwrite: overwrite the MESS file currently in the file system
nprocs: number of processers to use during the MESS calculation
inpname:

Transport Driver Tasks
----------------------

'onedmin'

Process Driver Tasks
--------------------

'freqs'
'energy' 
'geo':
'zmatrix':
'enthalpy'
'coeffs': 


