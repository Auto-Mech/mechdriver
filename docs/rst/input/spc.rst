
species.dat
===========

    List of species to perform calculations for.

    These species may also correspond to those in a chemical mechanism provided in the input file: `mechanism.dat`.

    Input provided in the form of a comma-separated value (csv) file.

    General format::
        <header1>,<header2>,...,<headerN>
        <spc1-inf1>,<spc1-inf2>,...,<spc1-infN>
        <spcN-inf1>,<spcN-inf2>,...,<spc1-infN>

    Allowed:
        (1) name
        (2) inchi
        (3) smiles
        (4) mult  [for multiplicity]
        (5) charge
        (6) sens. [for sensitivity]

    Note::
        Certain input (e.g.,InChI strings and SMILES strings) may have commas and other symbols which break the CSV file parser. To avoid this issue, these strings may be given in single quotes (‘).

    Certain information is requred in the input: name, mult, and ich/smil so that we can
    construct basic description of species

    Note that internally, all the species are represented by InChI strings. If only SMILES are provided for a given species, that SMILES string will be converted to an InChI string by the code. This is because these representations allow greater precision in structural description, such as for stereochemistry.

    Note::
        Certain input (e.g.,InChI strings and SMILES strings) may have commas and other symbols which break the CSV file parser. To avoid this issue, these strings may be given in single quotes (‘).


Species.dat
-----------

Auxiliary input for species that specifies where more complex info

Explanation of ts nomenclature

General format::

    spc
        key1=val1
        key2=val2...
    end spc

---somewhere explain lists

Names must be chemkin.csv, ts_name, global.

See global section.


.. csv-table:: keywords for spc
    :file: tables/spc_keys.csv
    :header-rows: 1
    :widths: 10, 10, 10, 10

Explantions:

(seperate into core, task specific, driver specific)

inchi: InChI string (up to stereochemical layer [ref])
smiles: SMILES string (any form)

mult: spin multiplicity

charge: electric charge

inchikey: InChI key

sensitivity: sensitivity value for mechanism

tors_names: list of dihedral coords that define internal rotors

sym_factor: overall symmetry factor (external+internal)

hind_inc: increment of rotor potential (in degrees)

mc_nsamp: number of samples for Monte Carlo sampling of conformers [give formula]

tau_nsamp: number of samples for Monte Carlo sampling of conformers [give formula]

etrans_nsamp: samples for OneDMin

bath: bath gas molecule species interacts with OneDMin

smin: minimum intermolecular distance for OneDMin

smax: maximum intermolecular distance for OneDMin

lj: lennard-jones params [used for mess writing]

edown: alpha, n for energy down model (mess writing)

active: active space variables 

zma_idx: index for zma in save filesystem

rxn_dirn: direction to look for TS in 

kt_pst: k(T) for PST theory

temp_pst: temperature for k(T)

n_pst: n parameter for PST model potential

pst_params: [pre-exponential factor, pre-exponential power] Assuming a PST model potential

Maybe include  examples for input?

