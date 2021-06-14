
Species
=======

Species information is provided in two files. One of which is required.


species.csv
-----------

    Enumerates the names of all mechanism species, alongside including basic necessary structural descripters. 

    Currently, the input correspond to the general format of a comma-separated value (csv) file::
        <header1>,<header2>,...,<headerN>
        <spc1-inf1>,<spc1-inf2>,...,<spc1-infN>
        <spcN-inf1>,<spcN-inf2>,...,<spc1-infN>

    Example (with required info)::
        name,inchi,mult,charge

    Note that certain pieces of information are futher enclosed in single quotes('), to not break the sections, particularly since they have commas.

    Internally, all mechanism species are represented with their given names and their InChI string representations. Therefore, these must be provided, or SMILES may also be provided as an alternative. The SMILES will be converted interanlly.

    Same Example with SMILES::
        name,inchi,mult,charge

    Pre-procesing can be done to add sterochem inchi and thermo species (see section).


Species.dat
-----------

Auxiliary input for species that specifies where more complex info for species as well as potential transition states.

The transition state name::
    ts_<pesidx>_<channeidx>_<tsidx>

where <pesidx> and <channelidx> correspond to the indices for the reaction PESs and channels defined in the mechanism file. 

The <tsidx> is the filesystem idx used to different different structures for a given transition state.

General format::

    spc <spc_name>
        key1=val1
        key2=val2...
    end spc


.. csv-table:: keywords for spc
    :file: tables/spc_keys.csv
    :header-rows: 1
    :widths: 10, 10, 10, 10


..
seperate or describe keys based on core definitions and ones that are task specific and driver specific. See global section.

