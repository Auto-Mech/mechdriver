
run.dat
-------

Central input file which describes all of the individual calculations the user wishes to run for one instance of MechDriver. User-defines these calculations via several containing lists of various formats that are detailed below.  Information in these sections will rely on definitions the user provides in the other input files.  These sections fall under three broad categories:

|    (1) Options for various I/O :ref:`input<input-block>`
|
|    (2) Lists of Reactions and Species defined in mechanism to run calculations for
|            pes, spc
|
|    (3) Lists of Tasks for Driver Tasks (electronic structure, kinetics, thermo methods):
|            els, thermo, trans, ktp, process

As can be seen there are multiple possible sections the user may define from type (2) and (3). Each of these sections are defined broadly via::

    <block_name>
        <block_contents>
    end <block_name>


.. _input-block:

Input Section
~~~~~~~~~~~~~

Specifies various I/O options such as the paths where calculations are performed and saved, as well as what the input mechanism an species file types are.

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


Chemistry Sections
~~~~~~~~~~~~~~~~~~


Tells the code what to do calculatiions to do for certain parts of the mechanism defined in the mech.dat file and spc.da file.

Currently both can be given in an output; however the PES block is given by
priority and the SPC block is combined with it. Certain cases where a block
is required.

All of the components listed in the two blocks will be done the instance of mechdriver.

Lists the Potential Energy Surfaces and constituent channels
(composed of species and transition states) for which to perform calculations.

The indices for the PESs and channels correspond to the those provided in input
`mechanism.dat` file. Similar to the PES Section, this section lists the species for which to perform calculations.
The indices for the species correspond to the order of the species provided in the input `species.csv` 

Note that multiple PES-Channel and species combinations can be given in one input. Lists can be extended across multiple lines. All of the components listed in the two blocks will be done the instance of mechdriver.

Note that all drivers can take both PES and spc lists. Will do them in that order. However, kTPDriver requires PES lists to be given.

General Format::

    spc
        <Species indices>
    end

    pes
        <PES indices>: <Channel indices>
    end

Example 1::

    pes
        1: 1          # Runs Channel 1 on PES 1
        2-4: 1-2      # Runs Channels 1,2 on PES 2,3,4
        5-8: 1,6-10   # Runs Channels 1,6,7,8,9,10 on PES 5,6,7,8
        9,12: 1-4,6   # Runs Channels 1,2,3,4,6 on PES 9,12
    end

Example 2::

    spc
        1          # Runs species 1
        1,3        # Runs species 1,3
        1-4        # Runs species 1,2,3,4
        1,3-5      # Runs species 1,3,4,5
        1-3,6      # Runs species 1,2,3,6
        1-3,5-7    # Runs species 1,2,3,5,6,7
    end


Task Sections
~~~~~~~~~~~

For every task block defined in the input, this will signal MechDriver to run its various sub-Drivers.

General format::

    tsks
        <object>  <task>  <keyword1=value  keyword2=value …>
        <object>  <task>  <keyword1=value  keyword2=value …>
        <object>  <task>  <keyword1=value  keyword2=value …>
        …
    end

| where 
| 
| <object> is either spc, ts. Optional for certain drivers
| <task> is what electronic structure calculation to be run on object.
keyword=value cannot have spaces in between them.

Each task is given in the following format <obj>_<job>


Comments
~~~~~~~~

    In general, all text preceded by `#` symbols will be ignored bythe parser. 

    As a trick for commenting out entire sections, comment out the header line of the section that will cause the entire section to be ignored during parsing.

This is an easy way to turn off an entire driver without commenting out several lines.

Commenting out entire sections::

    # section
        dwdqwd
        fwefdv
        fwefwe
    end

