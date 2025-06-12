AutoIO
=========

*A Package of the AutoMech Suite*

Andreas V. Copan, Kevin B. Moore III, Sarah N. Elliott, and Stephen J. Klippenstein

--------------------------------------------------------------------------------------

.. toctree::
    :glob:
    :maxdepth: 2

    index
--------------------------------------------------------------------------------------

Overview
~~~~~~~~

An essential feature of the  AutoMech suite is its seamless workflow through many types of calculations
including electronic structure, transport, master-equation, and modeling.
The AutoIO library manages the input and output operations and the interfaces to non-python codes
to make this workflow possible.
The AutoIO library contains two conda packages: autoio-base and autoio-interfaces. AutoIO-base is the lowest-level
package in the AutoMech suite, whereas AutoIO interfaces in dependent on AutoChem.

Base
^^^^

The base input/output operations for AutoMech are:

* autoparse: a set of wrappers for python.re that allow for more intelligible regular expressions syntax.

* autoread: reads files and converts strings into AutoMech obects

* autowrite: converts AutoMech objects into strings and writes them

* ioformat: formats strings

Interfaces
^^^^^^^^^^^

We have built interfaces to the following codes:

.. list-table:: 
   :widths: 10 11
   :header-rows: 1

   * - Internal Codes
     - External Codes
   * - MESS
     - chemkin
   * - onedmin 
     - pac99
   * - projrot?
     - polyrate
   * - intder 
     - 
   * - thermp
     -
   * - varecof
     -

The AutoIO-Interfaces package also houses the autorun module which, after creating input for codes listed above, proceeds to execute
the requested task and parse its output.

A highlight of the AutoIO-interfaces package is elstruct module which, itself, has many interfaces to electronic structure codes.  
It can automatedly run the following calculations on the following external quantum chemistry softwares:

* Single point calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
* Geometry calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
* Gradients calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
* Frequency calculation: CFOUR(v2), Gaussian(v09, v16), Molpro(v2015), Orca(v4), PSI4
* VPT2 anharmonic analysis: Gaussian(v09, v16)
* IRC: Gaussian(v09, v16), PSI4

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

    >>> conda install autoio-base -c auto-mech 
    >>> conda install autoio-interfaces -c auto-mech

For users new to conda, we have :ref:`conda-instructions`.
Each AutoMech package is also available on `GitHub`_.

.. _GitHub: https://github.com/Auto-Mech/autoio



Tutorial
^^^^^^^^
The first step is to make sure the installation was successful by importing some of the modules in each package

.. code-block:: python

    >>> import autoread
    >>> import elstruct


Then we can move on to using the autoio modules:

* Base Tutorial\: :ref:`base-tutorial-hub`
    * :ref:`ioformat-tutorial-doc`
    * :ref:`autoparse-tutorial-doc`
    * :ref:`autoread-tutorial-doc`
    * :ref:`autowrite-tutorial-doc`

* Interfaces Tutorial\: :ref:`interfaces-tutorial-hub`
    * :ref:`ioformat-tutorial-doc`
    * :ref:`chemkin-tutorial-doc`

Documentation
~~~~~~~~~~~~~
    .. toctree::
        :maxdepth: 3
    
        submodule_base
        submodule_interfaces
