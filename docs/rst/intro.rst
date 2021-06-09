
Overview
========

MechDriver serves as the automated workflow code of the AutoMech suite. 

Its primarily utility is to facilitate the calculation of thermochemical, kinetic, and transport
parameters required to construct reliabel reaction mechanism models.

The design principle of MechDriver is such that to act as a workflow manager which interfaces to any number of Python and other compiled codes to calculate quantities of interest.


The user supplies (1) a mechanism containing reaction channels and species and
(2) desired property of interest.


MechDriver will take these input options and schedule a set of tasks to calculate these quantities.

In a simplified manner, it 


Citation
========

The following citations introduce and utilize MechDriver:

.. bibliography:: refs.bib
    :list: enumerated
    
    automech1
    automech2


Support
=======

Technical Details and Installation:
Kevin Moore [kmoore@anl.gov]

