
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


Chemistry Sections
~~~~~~~~~~~~~~~~~~

Two separate sections, a `pes` block and `spc` block can be placed in the input to 
detail the potential energy surfaces and species from the input mechanism (defined in the mechanism.dat and species.csv) for which to perform calculations.

General Format::

    pes
        <PES indices>: <Channel indices>
    end

    spc
        <Species indices>
    end

The indices for the PESs and channels correspond to the those provided in input
`mechanism.dat` file. Similar to the PES Section, this section lists the species for which to perform calculations.
The indices for the species correspond to the order of the species provided in the input `species.csv` 

PES Block Example::

    pes
        1: 1          # Runs Channel 1 on PES 1
        2-4: 1-2      # Runs Channels 1,2 on PES 2,3,4
        5-8: 1,6-10   # Runs Channels 1,6,7,8,9,10 on PES 5,6,7,8
        9,12: 1-4,6   # Runs Channels 1,2,3,4,6 on PES 9,12
    end

SPC Block Example::

    spc
        1             # Runs species 1
        2,5           # Runs species 1,2
        7-9           # Runs species 7,8,9
        10,13-15      # Runs species 10,13,14,15
        17-19,22      # Runs species 17,18,19,22
        24-26,28-30   # Runs species 24,25,26,28,29,30
    end

Note that any combination of PESs and channels can be given, as well as any number
of species. Moreover the user can define sets of indices on as many lines as they see fit. In the above example, a requested driver will exectute requested tasks for ALL the
channels and species detailed in the accompanying comment lines. Any redundancy of
among the blocks may cost the user time but generally will not breakt things because
of the interaction with the run-save filesystem.

The requirements for `which` block must be provided depends somewhat on the driver
the user wishes to use in the MechDriver run; however most drivers can utilize both
sections with:
ESDriver: uses either/both pes and spc
ThermoDriver: uses either/both pes and spc (ignores TS)
kTPDriver: uses only pes
TransDriver: uses either/both pes and spc
ProcDriver: ???

Again from the above, we see that most drivers can utilize one or the other blocks.
Meaning you can run a set of tasks having only defined one of the blocks. In fact,
most blocks can loop over species specified in both blokcs. In the cases where user
has provided both a pes and spc block for runs a driver that uses both. The driver
will first loop over all tasks from the PES block first, then loop over all tasks
for the SPC block.


Task Sections
~~~~~~~~~~~~~

For each subdriver the user wishes to run with MechDriver, the user provides a section
with a list of tasks they wish the driver to execute. If no subdriver section is provided,
the driver will not be run during the MechDriver instance. The format of the section
varies upon which driver is specified, although there is a general format::

    <driver name>
        <object>  <task>  <keyword1=value  keyword2=value …>
        <object>  <task>  <keyword1=value  keyword2=value …>
        <object>  <task>  <keyword1=value  keyword2=value …>
        …
    end <driver name>

The initial <object> is only utilized by ESDriver. Here the user specificies if the
user wishes to run a given task for `spc` or `ts` or `all`. For all drivers only
the task and keyword-value pairs are required.

Each task and keyword-value pair are driver specific. Moreover, only keywords are only
pertinent to certain tasks. See link? The only formatting
guidelines of note is that the the keyword-value pairs must be input as presented, where
each pair has only an "=" sign and no whitespace.

Based on the task, certain keywords may be required, but in general, the user may define
any number of allowed task keywords on the line.

Example::
    <es>

Note the names on `lvl` and `mod` keywords. These correspond to names defined in
theory.dat and models.dat.

Similar to the other blocks, any number of tasks can be provided be provided in the
subdriver task block.

While the subdriver task blocks can placed in the input file, the subdrivers are
launched in a specific order:
(1) ESDriver, (2) ThermoDriver, (3) TransDriver, (4) kTPDriver, (5) ProcDriver


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

