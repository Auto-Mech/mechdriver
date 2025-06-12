.. _preproc:


Pre-Process Mechanism
---------------------

MechDriver will, in principle, be able to perform with any mechanism. However, 


Files
~~~~~

    mechanism.dat
    species.csv
    sort.inp (maybe)


Functionality:

(1) Inclusion of Stereochemistry. Generates InChI strings for all species with full
    stereochemical specification. These are used for SPC and RXN filesystem storage.

Produce species.csv

(2) Adding thermochemical basis functions: Relative energies used thermochem and rate
calculations are determined using CBH-n style schemes requiring basis functions. Electronic
structure calculations must be performed for these as well, hence there is utility in having
these added to species.csv to be placed in task lists.

(3 Sorting/Filter Reactions:
 

Post-Process Mechanism
----------------------

There are several useful analytical processes that can be conducted on mechanism files
that are either obtained externally or built by MechDriver.

Plotting
~~~~~~~~

We have developed tools for reading the thermochemical and kinetic parameters
from a mechanism file and plot corresponding values at certain temperatures and pressures.

<img of plot>

As a part of this functionality, we are able to create comparison plots of species thermochem
and reaction kinetic values from two different mechanisms. The tools can use any mechanism, however,
we frequently use it to compare input mechanism to one built by MechDriver.

<img of comparison plot>


Checking
~~~~~~~~

Assesses whether the mechanism satisfies certain physical constraints

    Check the mechanism

