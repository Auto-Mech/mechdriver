.. _execution:

Running
========================

Execution of the code is very simple; the user only needs to run the mechdriver executable with python. The input files are a bit more complex, but lets make sure the executable runs first!

.. code-block:: console
    
   mdriver -i <path/to/input_directory> -o <name of output file>

It should be emphasized that not only do all of the AutoMech dependencies need to be handled (by Conda or direct install), but any electronic structure package that will be used in the calculations

Depending on how dependencies are set-up, a simple BASH script might look at like

.. code-block:: bash

    #!/bin/bash
  
    # Set path
    export OUTPUTFILE=$2
    
    # SSH into node and run calculation
    ssh -n $1 " conda activate <AutoMech Env Name>                                          ;
                cd $PWD                                                                     ;
                module load <electronic structure package 1>                                ; 
                module load <electronic structure package 2>                                ; 
                mdriver -i $PWD -o $OUTPUTFILE &                                            ;
                disown %1 "

This would SSH onto a node, load the Conda environment to handle the MechDriver dependencies, load electronic structure codes, then execute MechDriver to pipe output to a file. The specific code will differ based on the user input (either proceeding through SSH to a compute node, or proceeding through a queiing system).

Currently, no multicore parallelilization is utilized directly by MechDriver, and processors for codes are defined in the input files. 

(Add argparse options to mechdriver)

As a Python code, individual modules and functions can be imported for use in other codes; however, this may not be the best use of this code. Most functionality of interest to import will likely exist in the lower-level libraries of the AutoMech suite.

.. note::
    Up next: 
    :ref:`execution`.

..

.. note::
   Return to:
   :ref:`manual`.

