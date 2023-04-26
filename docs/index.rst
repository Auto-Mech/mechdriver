.. MechDriver documentation master file

**********
MechDriver
**********

*A Package of the AutoMech Suite*

Sarah N. Elliott, Kevin B. Moore III, Andreas V. Copan 
Clayton Mulvihill, Luna Pratali Maffei, and Stephen J. Klippenstein


Overview
~~~~~~~~

MechDriver is the capstone of the AutoMech suite. It drives elaborate workflows of modules from
the AutoMech package to produce electronic structure, thermochemical, kinetic, and transport properties from
simple mechanism inputs.  These parameters are formatted into standard functional forms used
in the simulation of reaction mechanism models.  Follow these links for descriptions of the 
input files, a glossary of the keywords, and for example input files. The code structure and install
instructions are described below:

--------------------------------------------------------------------------------------

.. toctree::
    :maxdepth: 2
    :caption: Contents:

    rst/usersmanual


--------------------------------------------------------------------------------------

mechdrivers
^^^^^
The main automech.py code will drive any of the five drivers of the code, dependent on the user's input,
on a provided mechanism.  These drivers are:

* ESDriver: carries out electronic structure tasks for requested species and theories
* ThermoDriver: carries out partition function calculations and tranformations for requested species and models
* kTPDriver: carries out T,P-dependent rate calculations and fitting for multichannel potential energy surfaces 
* TransDriver: carries out energy transfer calculations for a given mechanism
* ProcDriver: post-processing module to format data

mechroutines
^^^^
To fascilitate the drivers, we have created subsets of routines

* es: modules to set up and run electronic structure jobs like geometry optimization, frequency analaysis, conformer searching, angle scanning, single point calculations, 
* ktp: organizes reaction, species, and transition state info to construct indepndent PESes 
* models: gathers the partition function parameters for each stationary point in the mechanism, and applies any user requested scaling routines
* proc:  querries the database for user requested properties
* thermo: produces the partition functions and transfroms it into thermchemical properties and NASA polynomials 
* trans:


mechlib
^^^^^^^^^^^
Library to house small, specific functions used by both the overrall drivers and their subroutines.
* amech_io: parsing functions specific to mechdriver user input and writing functions for code output

* filesys: code for reading and saving specific molecular information to the filesystem, and well as some basic transformations for such data, and locator functions to determine where to read/write from

* reaction: sets up reaction information


Each of these subdrivers function by making appropriate calls to lower level Python libraries packaged
in the AutoMech suite:

    | `autochem: <https://sne-autochem.readthedocs.io/en/latest/>`_ molecular representation  
    | `autoio: <https://sne-autoio.readthedocs.io/en/latest/>`_  I/O parsing and execution for external codes
    | `autofile: <https://sne-autofile.readthedocs.io/en/latest/>`_ builds and manages the RUN/SAVE filesystem 
    | `mechanalyzer: <https://mechanalyzer2-kev.readthedocs.io/en/latest/>`_ pre- and post-processing mechanisms; handles mechanism objects internally

AutoMech is also distributed with binary codes:

    | MESS
    | ThermP+PAC99
    | OneDMin
    | PIPPy

MESS is particularly important and forms the backbone of the AutoMech code as it handles all of the 
master equation rate constant expressions.

Getting Started
~~~~~~~~~~~~~~~
Installation
^^^^^^^^^^^^^
We have conda packages on the anaconda cloud for all of our packages. To install them,
set up an environment for AutoMech.  You can use the environment we have prepared for the
suite auto-mech-env.  Then activate your environment and install the autoio packages.


.. code-block:: python

    >>> conda env create auto-mech/auto-mech-env
    >>> conda activate auto-mech-env

For users new to conda, we have :ref:`conda-instructions`.
Each AutoMech package is also available on `GitHub`_.

.. _GitHub: https://github.com/Auto-Mech/autoio



Tutorial
^^^^^^^^
* Base Tutorial\: :ref:`base-tutorial-hub`
    * :ref:`ioformat-tutorial-doc`


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Support
=======

| Sarah Elliott [elliott@anl.gov]

