
Code Scripts
------------

Pre:Sorting
~~~~~~~~~~~

Sorting: Able to sort the reactions according to various groupings.

Will describe the PES-Channel numbering system that is assigned internally and put in the output file. Helps
with running mechdriver.


.. code-block:: console
    sort_mech -i


mechanism.dat file::

    <STRING FOR MECHANISM>

species.csv file::

    <STRING FOR SPECIES file>

sort.dat file::

    isolate_submech # this section is not mandatory
        species = (
            C2H4
        )
    end
    
    sort_mech
        criteria = (
            subpes
            molecularity
            rxn_class_broad
        )
        n_criteria_headers = 0
    end

The isolate_submech section will <insert description>


Pre:Species Modification
~~~~~~~~~~~~~~~~~~~~~~~~

Much of the functionality requires all species to have assigned InChI strings with full stereochemistry.
Additionally, much of the energy calculations accesses CBH calculations which require basis species.

Including all of this information in an initial species.csv file is time-consuming and often does not
come from mechanisms provided externally. As such, we have 

.. code-block:: console

    sort_mech -i


Post:Rate Fitting
~~~~~~~~~~~~~~~~~

Able to fit all of the reactions that are given in a MESS output file to some desired functional form and output a CHEMKIN file.

The only required file for the script is a MESS output file from a rate calculation

One can also provide an auxiliary input file to convert the names from the MESS output (since it is common to use
simplified names for MESS files). This file is also


Note for JSON format, all of the strings for dictionary keys and values should be in double quotes.

One can also provide a file to convert the names of the species in the MESS file to that of the
Can also take a label file

.. code-block:: console

    sort_mech -i

Example::
    
    <MESSRATE output>
    <label eux file>

Will read all of the rates and filter them


Post:Plotting
~~~~~~~~~~~~~

Plot rate constants from one or more mechanisms.

Note unless plotting at few pressures, it is not recommended to plot more than two mechanisms.

required files:
mechanism.dat files
therm.dat files

.. maybe just have single mechanism 

