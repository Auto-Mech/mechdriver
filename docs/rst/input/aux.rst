
auxiliary
=========

geom.xyz
--------

For species where the user wishes to seed some mechanism process with an input geometry, the user can supply any number of geom.xyz
files. The format is a standard xyz file except the comment line means something

<natoms
<name of species in species.csv>
<xyz coordinates>

The key note is the comment line. The name used here must match some name of the species.csv file.

Can also be a transition state

this will be used


files.aspace
------------

For multireference calculations, the user may specify a wavefunction guess template for a molpro calculation. This will be inserted in front of wfn calculations.

example:

! <transition state name>
<rest of the wfn template>

Internally, we have only very mild abilities to create guess wavefunctions. We can make
(2,2) space assuming the two highest occupied orbitals are the two reactive radical orbitals.

