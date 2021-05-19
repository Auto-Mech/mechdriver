
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


.. list-table:: Keywords for the input section
   :widths: 25 15 25 50
   :header-rows: 1

   * - Keyword
     - Required
     - Allowed
     - Default
   * - run_prefix
     - x
     -
     - None
   * - save_prefix
     - x
     -
     - None
   * - inp_mech
     - x
     - chemkin
     - chemkin
   * - out_mech
     - x
     - chemkin
     - chemkin
   * - inp_spc
     - x
     - csv
     - csv
   * - out_spc
     - x
     - csv
     - csv


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

Allowed Tasks
------------

Each task is given in the following format <obj>_<job>
 
'init_geom'
'find_ts'
'conf_pucker'
'conf_samp'
'conf_energy'
'conf_grad'
'conf_hess'
'conf_vpt2'
'conf_prop'
'conf_opt
'hr_scan'
'hr_grad'
'hr_hess'
'hr_energy'
'hr_vpt2'
'hr_reopt': (',)),
'tau_samp': (
'tau_energy':
'tau_grad': (
'tau_hess': (),
'rpath_scan':
'rpath_energy,
'rpath_grad':
'rpath_hess':
# Transport Driver Tasks
'onedmin': (('spc',), (BASE + TRANS)),
# Process Driver Tasks
'freqs': (('spc', 'ts', 'vdw'), PRNT + ('scale',)),
'energy': (('spc',), PRNT),
'geo': (('spc',), PRNT),
'zmatrix': (('spc',), PRNT),
'enthalpy': (('spc',), PRNT),
'coeffs': (('spc',), ()),
# KTP/Therm
'write_mess': ((), ('kin_model', 'spc_model', 'overwrite')),
'run_mess': ((), ('kin_model', 'spc_model', 'nprocs', 'inpname')),
'run_fits': ((), ('kin_model',)),


Electronic Structure Driver Task Block
---------------------------------------

Specifies electronic structure tasks

Order matters


kTPDriver Task Block
--------------------

Order does not matter

