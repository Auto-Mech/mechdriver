Autofile
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

Autofile directs the filesystem framework for the AutoMech suite.  It uses a dual architecture: the run/save filesystems.
Under each, Autofile defines rigid Dataseries trunks, branches, and leaves that each store specific data. The save system is the
database for AutoMech, storing molecular properties like electronic energies, hessians, and torsional profiles as well as
and transition state information like transitory bonds and reaction types.

.. image:: autofile.png
  :width: 800

Getting Started
~~~~~~~~~~~~~~~
Installation
^^^^^^^^^^^^^

We have conda packages on the anaconda cloud for all of our packages. To install them,
set up an environment for AutoMech.  You can use the environment we have prepared for the
suite auto-mech-env.  Then activate your environment and install the autofile package.


.. code-block:: python

    >>> conda env create auto-mech/auto-mech-env
    >>> conda activate auto-mech-env

    >>> conda install autofile -c auto-mech

For users new to conda, we have :ref:`conda-instructions`.
Each AutoMech package is also available on `GitHub`_.

.. _GitHub: https://github.com/Auto-Mech/autofile

Tutorial
^^^^^^^^
The first step is to make sure the installation was successful by importing autofile

.. code-block:: python

    >>> import autofile


Then we can move on to using the autofile module:

* :ref:`spc-tutorial-doc`
* :ref:`thy-tutorial-doc`
* :ref:`cnf-tutorial-doc`
* :ref:`zmat-tutorial-doc`
* :ref:`scn-tutorial-doc`
* :ref:`rxn-tutorial-doc`
* :ref:`ts-tutorial-doc`


Documentation
~~~~~~~~~~~~~
    .. toctree::
        :maxdepth: 4

        submodule_fs
        submodule_model
        submodule_schema
        submodule_data_types
        submodule_info
        submodule_io
        submodule_json


