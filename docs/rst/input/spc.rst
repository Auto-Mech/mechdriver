
species.csv
-----------

Enumerates the names of all mechanism species, alongside including basic necessary structural descripters. 

Currently, the input correspond to the general format of a comma-separated value (csv) file::
    <header1>,<header2>,...,<headerN>
    <spc1-inf1>,<spc1-inf2>,...,<spc1-infN>
    <spcN-inf1>,<spcN-inf2>,...,<spc1-infN>

See <csv link> for common practices for CSV files

Example (with required info)::
    name,inchi,mult,charge

Note that certain pieces of information are futher enclosed in single quotes('), to not break the sections. InChI and 
SMILES commonly have characters (e.g., commas) that will break the CSV parser.

Internally, all mechanism species are represented with their given names and their InChI string representations. It is 
thus recommended to run MechDriver with a species.csv where all species have full stereochemically labeled InChI strings.
This can be done with pre-processing modules found in MechAnalyzer (see section). Similarly, all CBH species should be added.

However, since InChI strings are difficult to read, MechDriver can take SMILES as an alternative and convert then internally
during runtime. Once of the these two must be present alongside the name and mult.

Same Example with SMILES::
    name,inchi,mult,charge


species.dat
-----------

Auxiliary input for species that specifies where more complex info for species as well as potential transition states.
Here, various sections of keywords are provided for species and transition states that comprise the mechanism. These 
sections are each formatted like as::

    spc <spc_name>
        key1=val1
        key2=val2...
    end spc

where <spc_name> corresponds to the species the user wishes to set the keyword values to. For species, these names
must correspond to those defined in the species.csv file (and by extension the mechanism.dat file). For transition states,
we have devised a specific nomenclature used throughout the code::

    ts_<pesidx>_<channeidx>_<tsidx>

where <pesidx> and <channelidx> correspond to the indices for the reaction PESs and channels defined in the mechanism file.
The <tsidx> is the filesystem idx used to different different structures for a given transition state.

For example, the isomerization of 1-butyl radical to 2-butyl radical (CH3CH2CH2CH2 <-> CH3CH2CHCH3) can occur via two distinct
transition states, which will be given tsidx of 0 and 1. These labels are arbitrary. Generally, setting tsidx to 0 in the name
will be sufficient.

Of course, there may be situations where the user wishes to set a parameter to a non-default value for several hundereds or thousands of species. In this case specifying a species block for all of these and setting names would be too cumbersome. For these cases, a
user may define a `global` block::

    spc global
        key1=val1
        key2=val2...
    end spc

where the keyword-value parameters specified here will be applied to ALL species in the mechanism. Note that a global block
can be defined alongside any number of specifically named blocks , like so::

    spc global
        key1=False
    end spc
    
    spc CH4
        key1=True
    end spc

    spc C2H6
        key1=True
    end spc

in this case all x species in the mechanism will bhave key1 = False except CH4 and C2H6 where key1 is True. In essence,
named spc block parameters will override those given in the global block. This sets a hiearchy of keyword setting to be:
(1) params defined in named, specific species blocks
(2) params defined in the global block (if one exists)
(3) params defined by internal defaults.

