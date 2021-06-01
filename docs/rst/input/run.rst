
run.dat
=======

Central input file which details the calculations that will be run during an instance of MechDriver. The user supplies several aspects of the calculation with sections::

    <block_name>
        <block_contents>
    end <block_name>

Three different primary blocks set a full mechdriver run.

    (1) Input/Output:
    (2) Chemistry Lists (parts of mechanism to run):
            pes, spc
    (3) Driver Tasks (electronic structure, kinetics, thermo methods):
            els, thermo, trans, ktp, process

One or more sections from type (2) and (3) can be provided.

For many of the above, the user will supply labels that correspond to blocks of info that is pulled from the other input file. Auxlilary information and defintions will come from other input


Input Block
-----------

Simple block that specifies forms of input

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

Tells the code what to do calculatiions to do for certain parts of the mechanism defined in the mech.dat file and spc.da file.

Currently both can be given in an output; however the PES block is given by
priority and the SPC block is combined with it. Certain cases where a block
is required.

All of the components listed in the two blocks will be done the instance of mechdriver.

Lists the Potential Energy Surfaces and constituent channels
(composed of species and transition states) for which to perform calculations.

The indices for the PESs and channels correspond to the those provided in input
`mechanism.dat` file. Similar to the PES Block, this section lists the species for which to perform calculations.
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
        1: 1         # Runs Channel 1 on PES 1
        1-2: 1-2     # Runs Channels 1,2 on PESs 1,2
        1,3: 1-4     # Runs Channels 1,2,34 on PESs 1,3
        1,3: 1-4,6   # Runs Channels 1,2,34 on PESs 1,3
    end

Example 2::

    pes
       1: 1-5        # Runs Channels 1,2,3,4,5 on PES 1
       2: 1,6-10     # Runs Channels 1,6,7,8,9,19 on PES 2
       5: 8          # Runs Channel 8 on PES 5
    end

Example 3::

    spc
        1          # Runs species 1
        1,3        # Runs species 1,3
        1-4        # Runs species 1,2,3,4
        1,3-5      # Runs species 1,3,4,5
        1-3,6      # Runs species 1,2,3,6
        1-3,5-7    # Runs species 1,2,3,5,6,7
    end


Task Blocks
----------

For every task block defined in the input, this will signal MechDriver to run its various sub-Drivers.

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


Commenting out Blocks
---------------------

    In general, all text preceded by `#` symbols will be ignored bythe parser. 

    As a trick for commenting out entire sections, comment out the header line of the section that will cause the entire section to be ignored during parsing.

This is an easy way to turn off an entire driver without commenting out several lines.

Commenting out entire sections::

    # section
        dwdqwd
        fwefdv
        fwefwe
    end

