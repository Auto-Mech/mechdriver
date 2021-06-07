
Auxiliary Input Files
=====================

XYZ Geometry Files
------------------

For species where the user wishes to seed some mechanism process with an input geometry, the user can supply any number of geom.xyz
files. The format is a standard xyz file except the comment line means something

File format::

<number of atoms>
<name of species in species.csv file>
<xyz coordinates>

The key note is the comment line. The name used here must match some name of the species.csv file.

Can also be a transition state

this will be used


Active Space Template Files
---------------------------

For multireference calculations, This will be inserted in front of wfn calculations.

File format::

    ! <transition state name>
    <rest of the wfn template>

Internally, we have only very mild abilities to create guess wavefunctions. We can make
(2,2) space assuming the two highest occupied orbitals are the two reactive radical orbitals.

Currently, only Molpro is supported.

