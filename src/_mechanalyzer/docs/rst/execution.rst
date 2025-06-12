
Installation and Running
========================

Installation
------------

Conda

Activate environment or install mechanalyzer and dependencies directly

.. code-block:: console

    conda install -c auto-mech mechanalyzer

This will install mechanalyzer as well as the accompanying submodules ratefit and thermfit.


Running
-------

Most immediately useful functionality of MechAnalyzer can be executed by standalone scripts existing in the
bin directories of the sub-modules the code. If in the Conda environment, these will exist in the user's $PATH.

Each of the scripts can be executed at the command line (assuming script is in the path)

.. code-block:: console

    <script-name> <options>

 To observe the available options for each script:

.. code-block:: console

    <script-name> --help

The base input required is generally a mechanism and/or species file. However, for scripts with more involved 
options, additional auxiliary input files will be required for proper execution.

.. include:: code_scripts.rst


Importing
---------

Of course all of the modules and functions code can be imported for use as a library.

.. code-block:: python

    >>> import mechanalyzer
    >>> import thermfit
    >>> import ratefit

