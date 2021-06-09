
Object Description
------------------

Following are standard descriptions of complex objects that are used throughout the workflow

spc_dct
Description: Dictionary of necessary structural and jobtime information for various species and transition states of the mechanism.
See x for description of allowed keywords used to build the spc dct

.. code-block:: python

    {
        spc_name: {
            inf: inf
            inf: inf},
        spc_name: {
            inf: inf
            inf: inf},
        ts_name: {
            inf: inf
            inf: inf},
    }

thy_dct
Description: Dictionary of user-defined electronic structure methods that fully define methodology and jobtime run options.

.. code-block:: python

    {
        thy_name: {
            inf: inf
            inf: inf},
        thy_name: {
            inf: inf
            inf: inf}
    }

model_dct
Description: Dictionary of user-defined model parameters used to calculate the partition functions for thermochemistry master equation models for kinetics.

.. code-block:: python

    {
        model_name: {
            inf: inf
            inf: inf},
        model_name: {
            inf: inf
            inf: inf}
    }

pes_dct:
Description: Reaction Channels from mechanism file sorted into potential energy surfaces

.. code-block:: python

    {(pes_formula, pes_idx, subpes_idx): (
        (chnl_idx, (chnl_reacs,), (chnl_prods)), 
        (chnl_idx, (chnl_reacs,), (chnl_prods)), ...}

pes_rlst
Description: Dictionary of PESs to loop over in various Drivers. Mirror structure of pes_dct

.. code-block:: python

    {(pes_formula, pes_idx, subpes_idx): (
        (chnl_idx, (chnl_reacs,), (chnl_prods)), 

spc_rlst
Description: Dictionary containing species to loop over in various Drivers. Designed to mimic pes_rlst for code simplicity.

.. code-block:: python

    {('SPC', 0, 0): (spc_idx, spc_idx2, ...)} 

run_rlst
Description: Combination of the pes and spc rlsts

spc_queue
Description: List of species to loop over for a set of driver tasks.

.. code-block:: python

    (spc_name1, spc_name2, ...)

