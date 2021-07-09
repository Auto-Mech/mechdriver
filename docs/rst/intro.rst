
Overview
========

MechDriver serves as the automated workflow code of the AutoMech suite. Acts as a manager
to execute a series of tasks involving the calculation of chemical data, including electronic structure
data as well as paramters describing thermochemistry, kinetic, and transport. These parameters are
formatted into standard functional forms used in the simulation of reaction mechanism models. 

To accomplish each of the specified tasks, MechDriver calls the appropriate, constituent subdriver:
    (1) ESDriver: electronic structure calculations (geometries, frequencies, energies, etc)
    (2) ThermoDriver: species thermochemistry (enthalpies, entropies, etc)
    (3) kTPDriver: T,P-dependent rate constants of multichannel potential energy surfaces
    (4) TransDriver: energy transfer and chemical transport
    (5) ProcDriver: post-processing module to format data

Each of these subdrivers function by making appropriate calls to lower level Python libraries packaged
in the AutoMech suite:

    | `autochem: <https://sne-autochem.readthedocs.io/en/latest/>`_ molecular representation  
    | `autoio: <https://sne-autoio.readthedocs.io/en/latest/>`_  I/O parsing and execution for external codes
    | `autofile: <https://sne-autofile.readthedocs.io/en/latest/>`_ builds and manages the RUN/SAVE filesystem 
    | `mechanalyzer: <https://mechanalyzer2-kev.readthedocs.io/en/latest/>`_ pre- and post-processing mechanisms; handles mechanism objects internally

These libraries handle much of the complicated procedures for I/O of external codes and processing
data. MechDriver largely handles the sequence of tasks.

AutoMech is also distributed with binary codes:

    | MESS
    | ThermP+PAC99
    | OneDMin
    | PIPPy

MESS is particularly important and forms the backbone of the AutoMech code as it handles all of the 
master equation rate constant expressions.

Support
=======

| Kevin Moore [kmoore@anl.gov]
| Sarah Elliott [elliott@anl.gov]

