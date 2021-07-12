
Object Description
------------------

Following are standard descriptions of complex objects that are used throughout the workflow

Note that some objects can be constructed from objects of variable types. These are labeled as `obj`.


**species information dictionary (spc_dct)**
Description: Dictionary of necessary structural and jobtime information for various species and transition states of the mechanism.
See x for description of allowed keywords used to build the spc dct

.. code-block:: python

    # Description
    {
        spc_name_i: spc_dct_i,
        spc_name_j: spc_dct_j,
        ts_name_i: spc_dct_k
    }

    # Type
    dict[
        str: dict[str]
    ]

    # Example (see spc_dct_i and ts_dct for example of those
    {
        'C3H7(1)': spc_dct_i,
        'C3H7(2)': spc_dct_i,
        'ts_1_1_0': ts_dct
    }

spc_dct_i
Description: Dictionary of info for single species (sub-dict of spc_dct)

.. code-block:: python
    
    # Description
    {
        param_key_i: param_value_i,
        param_key_j: param_value_j,
    }

    # Type
    dict[
        str: dict[str: obj]
    ]

    # Example
    {
        'smiles': 'CC[CH2]'                     # basic info from .csv file
        'mc_nsamp': [True, [10, 3, 2, 1, 10]]   # extra runtime info needed by ESDriver 
    }

ts_dct
Description: variant of spc_dct_i for TS with special information required for that

.. code-block:: python

    # Description
    {
        param_key_i: param_value_i,
        param_key_j: param_value_j,
    }

    # Type
    dict[
        str: dict[str: obj]
    ]

    # Example
    {
        'zrxn': automol.reac.Reaction object
    }


thy_dct
Description: Dictionary of user-defined electronic structure methods that fully define methodology and jobtime run options.

.. code-block:: python

    # Description
    {
        thy_name_i: {
            param_key_i: param_value_i,
            param_key_j: param_value_j
        }
        thy_name_j: {
            param_key_i: param_value_i,
            param_key_j: param_value_j
        }
    }

    # Type
    dict[
        str: dict[str: obj]
    ]

    # Example
    {
        'program': 'molpro2015',
        'method': 'ccsd(t)-f12',
        'basis': 'cc-pvdz-f12',
        'orb_res': 'RR'
    }

pes_model_dct/spc_mod_dct
Description: Dictionary of user-defined model parameters used to calculate the partition functions for thermochemistry master equation models for kinetics.

.. code-block:: python

    # Description
    {
        model_name_i: model_dct_i,
        model_name_j: model_dct_j,
    }


    # Example
    

pes_model_dct_i
Description: information for a pes model containing information needed to run and fit kinetics and thermo

.. code-block:: python

    {
        'pressures': pressure lst
        'temps': pressure lst
        
    }

spc_model_dct_i
Description: information on how to compute the partition functions

.. code-block:: python

    {
        'pf_comp_i' = {
            'mod': model_name
            'lvl': theory_namr
        }
        'pf_comp_j' = {
            'mod': model_name
            'lvl': theory_namr
        }
    }


pes_dct:
Description: Reaction Channels from mechanism file sorted into potential energy surfaces

.. code-block:: python

    # General description
    {
        (pes_formula, pes_idx, subpes_idx): (
            (chnl_idx, (chnl_reacs,), (chnl_prods)), 
            (chnl_idx, (chnl_reacs,), (chnl_prods)), ...
    }

    # Types
    dict[
        (str, int, int): (
            (int, tuple(str), tuple(str)),
            (int, tuple(str), tuple(str)),
        )
    ]

    # Example
    {
        ('C3H7', 0, 0): (
            (0, (('C3H7(1)',), ('C3H7(2)',)),
            (1, (('C3H7(1)',), ('C3H6', 'H')),
            (2, (('C3H7(2)',), ('C3H6', 'H')),
        )
        ('C3H9O', 1, 0): (
            (0, (('C3H8', 'OH'), ('C3H7(1)', 'H2O')),
        )
        ('C3H9O', 1, 1): (
            (0, (('C3H8', 'OH'), ('C3H7(2)', 'H2O')),
        )
    }


pes_rlst
Description: Dictionary of PESs to loop over in various Drivers. Mirror structure of pes_dct

.. code-block:: python

    # Description
    {(pes_formula, pes_idx, subpes_idx): (
        (chnl_idx, (chnl_reacs,), (chnl_prods)), 

    # Type
    dict[(str, int, int): (
        (int, tuple(str), tuple(str))
    ]

    # Example

spc_rlst
Description: Dictionary containing species to loop over in various Drivers. Designed to mimic pes_rlst for code simplicity.

.. code-block:: python

    # Description
    {('SPC', 0, 0): (spc_name_i, spc_name_j, ...)} 

    # Type
    dict[('SPC', 0, 0): tuple(str)]

    # Example
    {('SPC', 0, 0): ('C3H8', 'C3H7(1)', 'OH', ...)} 

run_rlst
Description: Combination of the pes and spc rlsts

spc_queue
Description: List of species to loop over for a set of driver tasks.

.. code-block:: python

    # Description
    (spc_name_i, spc_name_j, ...)

    # Type
    tuple(str)

    # Example
    ('C3H8', 'C3H7(1)', 'OH', ...)

pf_filesys
Description: Contains autofile objects describing the locations and paths for where to read electronic
structure data from the SAVE filesystem that is processed to produce final data for building some portion
of partition function or related information:

.. code-block:: python

    {model_component_i: (autofile.locs, autofile.min, ..)  get rest


spc_info
Description: Bundle of basic information to describe the physical and electronic structure of a species. Can
be used to access the RUN/SAVE fileystem layer for species using autofile code. Enough info to generate a new geometry as well.

.. code-block:: python

    # Description
    (InChI string, Charge, Multiplicity)

    # Type
    tuple(str, int , int)

    # Example
    ('InChI=1S/H2O/...', 0, 1)

thy_info
Description: Bundle of basic information to describe the . Parts 2,3,4 can
be used to access the RUN/SAVE fileystem layer for theory using autofile code.

.. code-block:: python

    # Description
    (program, method, basis, combined orb ref label)

    # Type
    tuple(str, str , str, str)

    # Example
    ('molpro2015','ccsd(t)-f12', 'cc-pvdz-f12', 'RR')

mod_thy_info
Description: thy info object except the orbital label corresponds to appropriate label for program and species
multiplcity

.. code-block:: python

    # Description
    (program, method, basis, combined orb ref label)

    # Type
    tuple(str, str , str, str)

    # Example
    ('molpro2015','ccsd(t)-f12', 'cc-pvdz-f12', 'R')

rxn_info
Description:



Test glossary object:

.. glossary::

    Sphinx
      Sphinx is a tool that makes it easy to create intelligent and beautiful documentation. It was originally created for the Python documentation, and it has excellent facilities for the documentation of software projects in a range of languages.

    RST
      |RST| is an easy-to-read, what-you-see-is-what-you-get plain text markup syntax and parser system. It is useful for in-line program documentation (such as Python docstrings), for quickly creating simple web pages, and for standalone documents. |RST| is designed for extensibility for specific application domains. The |RST| parser is a component of Docutils.

    Sublime Text
      Sublime Text is a sophisticated text editor for code, markup and prose. You'll love the slick user interface, extraordinary features and amazing performance.

