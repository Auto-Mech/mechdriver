.. _input:

********************
Input File Structure
********************

Required Files
==============

    run.dat
    models.dat
    theory.dat
    species.csv
    mechanism.dat (for mechanisms)

Auxiliary Input Files:
    species.dat
    _.xyz
    _.aspace

Run.dat
=======

The central input where the user specifies what will be run by the program
	(1) objects, e.g. species, reactions, PESs
	(2) drivers and tasks
	(3) levels of theoretical treatment.

For many of the above, the user will supply labels that correspond to blocks of info that is pulled from the other input file.

The run.dat file is comprised of four sections:

	(1) input
	(2) pes/spc
	(3) els
 	(4) trans
    (5) thermo
    (6) ktp
    (7) process

where each section is specified in run.dat as

	section
		section input
	end

Input Block
-----------

Simple block that specifies

Example::
    input
	    mech = chemkin
	    run_prefix = <path/to/run/prefix>
	    save_prefix = <path/to/save/prefix>
    end inpiut

.. list-table:: Keywords for the input section
   :widths: 25 15 25 50
   :header-rows: 1

   * - Keyword
     - Required
     - Allowed
     - Default
   * - run_prefix
     - x
     -
     - None
   * - save_prefix
     - x
     -
     - None
   * - inp_mech
     - x
     - chemkin
     - chemkin
   * - out_mech
     - x
     - chemkin
     - chemkin
   * - inp_spc
     - x
     - csv
     - csv
   * - out_spc
     - x
     - csv
     - csv


Chemistry Blocks
----------------

Tells the code what to do calculations for.

Currently both can be given in an output; however the PES block is given by
priority and the SPC block is combined with it. Certain cases where a block
is required.

Example (put one with both PES and SPC)


Potential Energy Surface Block
------------------------------

Lists the Potential Energy Surfaces and constituent channels
(composed of species and transition states) for which to perform calculations.

The indices for the PESs and channels correspond to the those provided in input
`mechanism.dat` file. Note that multiple PES-Channel combinations can be
given in one input.

Required for:
    kTPDriver

General Format::
    pes
        <PES indices>: <Channel indices>
    end

Example 1::
    pes
        1: 1         # Runs Channel 1 on PES 1
        1-2: 1-2     # Runs Channels 1,2 on PESs 1,2
        1,3: 1-4     # Runs Channels 1,2,34 on PESs 1,3
        1,3: 1-4,6   # Runs Channels 1,2,34 on PESs 1,3
    end

Example 2::
    pes
       1: 1-5
       2: 1,2-10
       5: 8
    end

Species Block
------------

Similar to the PES Block, this section lists the species for which to perform calculations.

The indices for the species correspond to the order of the species provided in the
input `species.csv` file. Note that species can be extended across multiple lines

Required for:
    ThermoDriver

General Format::
    spc
        <Species indices>
    end

Example::
    spc
        1          # Runs species 1
        1,3        # Runs species 1,3
        1-4        # Runs species 1,2,3,4
        1,3-5      # Runs species 1,3,4,5
        1-3,6      # Runs species 1,2,3,6
        1-3,5-7    # Runs species 1,2,3,5,6,7
    end


Task Block
----------

The overall section is

tsks
	<object>  <task>  <keyword1=value  keyword2=value …>
    <object>  <task>  <keyword1=value  keyword2=value …>
    <object>  <task>  <keyword1=value  keyword2=value …>
	…
end

Above, the <object> is either spc (species, i.e. reacs, prods) or ts

Above the <task> is what electronic structure calculation to be run on object.

keyword=value cannot have spaces in between them.


Electronic Structure Driver Task Block
---------------------------------------

Specifies electronic structure tasks

Order matters


kTPDriver Task Block
--------------------

Order does not matter


Model.dat
=========

File used to describe how the partition functions will be constructed

kin model: specify the conditions by which rates./thermo is calculated and fit

    pressures: pressures used to calculate k(T,P) values with the master equation
    rate_temps: temperatures used to calculate k(T,P) values with the master equation
    therm_temps: temperatures used to calculate partition functions used for
        thermochemical parameter determination
    rate_fit: methodology used to fit the rate constants to functional forms
        fit_method: functional form to fit rate constants to
        pdep_temps: list of temperatures at which to assess pressure dependence
        pdep_tol: % error tolerance for determining if reaction is pressure dep
        pdep_pval: pressure value to get rate constants for if no pressure dependence
        pdep_plow: low pressure to assess pressure difference to
        pdep_phigh: high pressure to assess pressure difference
        arrfit_dbltol: % error tolerance for
        troefit_params: list of parameters to fit in Troe expression
    therm_fit: methodology used to fit thermochemical parameters to NASA polynomials
        ref_scheme: CBH scheme used to calculate energies
        ref_enes: reference energies used to calcualte delta H

spc model: specify the means for calculating the constituent components of the
           partition functions used in rate and thermochemistry calcualtion.

        ene: relative energies
        rot: rotational pf, geometry and
        vib: vibrational pf, frequencies
        tors: vibrational pf, internal rotor degrees of freedom
        symm: symmetry number
        ts: handle transition states

    mod: model used to calculate the
    geolvl: electronic structure level to read geometry
    enelvl: electronic structure level to read energy
    vpt2lvl: electronic structure level to get from VPT2 calc
    tunnel: tunneling model
    sadpt: Transition state theory for saddle points
    nobar: Transition state theory method used for barrierless reactions
    wells: manner by which build entrance-/exit-channel wells for a TS

    theory levels correspond to ones provided in the theory.dat file

for enes
Example::
    ene = [0.2, lvl_cc_df]



Theory.dat
==========

Details various electronic structure methods.

Theory levels are given as a series of blocks like so:

level <level_name>
		<keyword1 = value>
		<keyword2 = value>
  		…
end level

Above, <level_name> serves as a descriptor of the method that is used in other parts of the code. This string can be whatever the user wishes. It must contain no whitespace.

Any number of these level blocks can be given in theory.dat. The only blocks used by the code are those specified in run.dat and models.dat.


Mechanism.dat
=============

List of all the reactions that AutoMech may run.

Currently, the input file takes the form of a CHEMKIN input file. The full details of the format can be found in the CHEMKIN manual.
However, a simplified version will suffice for running AutoMech. The file can take the form::

    REACTIONS
      Reac1+Reac2=Prod1+Prod2         1.00 0.00 0.00
      Reac3=Prod3                     1.00 0.00 0.00
      Reac4=Prod3+Prod5               1.00 0.00 0.00
      …
    END

Above, are names of chemical species undergoing the reaction. These names must correspond to the names given in the species.csv file.

The numbers included in the above example are included for the input parser to work (but are also ignored). Note that rate-constant fit parameters normally included in the CHEMKIN file are ignored.

Reaction lines preceded by a ‘!’ are ignored by the parser. (maybe include #?)

The reactions do not need to be in any particular order, as AutoMech can sort the reactions into PESs itself.

The only reactions that will be run by the code, will be those specified in run.dat.


species.csv and species.dat
===========================

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

spc
    key1=val1
    key2=val2...
end spc

---somewhere explain lists

Names must be chemkin.csv, ts_name, global.

See global section.


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


global section
--------------

Many of the inputs can make use of a 'global' section that can set all keywo


********************
Full Minimal Example
********************

Input
=====

The following are a set of inputs for a minimal example treatment of the CH3OH+H reaction system. Here we use a relatively simple theoretical treatment:
    (1) electronic structure: low-level DFT methods
    (2) kinetics: fixed transition state theory with 1DHR torstion treatments

In this file, we have specified to calculate the geometries, vibrational frequencies, and 1-dimensional torsional potentials at the `lvl_wbs` level and the energies at the `cc_lvl_d`. These are defined in the theory.dat file.

To set up the chemical reactions and species for the input mechanism, we set

mechanism.dat file::

    REACTIONS
        <mech: propene+H?>
    END

species.csv file::

    name,smiles,mult,charge
       <mech spc>
       <include CH4 for thermo>

run.dat file::

    input
        inp_mech = chemkin
        inp_spc = csv
        run_prefix = /fake/path/to/run
        save_prefix = /fake/path/to/save
    end input

    pes
        1: x-y
    end pes

    spc
        <idx for CH4>
    end spc

    els
        spc init_geom     runlvl=wbs   inplvl=wbs
        ts  find_ts       runlvl=wbs   inplvl=wbs
        all hr_scan       runlvl=wbs   inplvl=wbs  tors_model=1dhrfa
        all conf_energy   runlvl=ccdz  inplvl=wbs
        all conf_hess     runlvl=wbs   inplvl=wbs
    end els

    thermo
        write_mess      kin_model=global  spc_model=global
        run_mess        kin_model=global  spc_model=global
        run_fits        kin_model=global  spc_model=global
    end thermo

    ktp
        write_mess      kin_model=global  spc_model=global
        run_mess
        run_fits        kin_model=global spc_model=global
    end ktp

Note that the pes specifies the global models. These models define the theoretical treatment used to build the MESS file and the rates

model.dat::
    kin global
        pressures = (
            0.1  1.0  10.0 100.0
        )
        rate_temps = (
            500. 600. 700. 800. 900. 1000.
            1100. 1200. 1300. 1400. 1500
            1600. 1700. 1800. 1900. 2000.
        )
        therm_temps = (
            200. 300. 400. 500. 600. 700. 800. 900. 1000. 1100. 1200.
            1300. 1400. 1500. 1600. 1700. 1800. 1900. 2000. 2100. 2200.
            2300. 2400. 2500. 2600. 2700. 2800. 2900. 3000.
        )
        rate_fit = (
            fit_method = plog
            pdep_temps = [500.0, 1000.0]
            pdep_tol = 20.0
            pdep_pval = 1.0
            arrfit_dbltol = 15.0
        )
        therm_fit = (
            ref_scheme = basic
            ref_enes = ANL0
        )
    end

    spc global
        ene = (
            lvl1 = ccdz
        )
        rot = (
            mod = rigid
        )
        vib = (
            mod = harm
            geolvl = wbs
        )
        tors = (
            mod = 1dhr
            enelvl = wbs
            geolvl = wbs
        )
        symm = (
            mod = sampling
            geolvl = wbs
        )
        ts = (
            tunnel = eckart
            sadpt = fixed
            wells = fake
        )
    end


In this file, we have specified to calculate the geometries, vibrational frequencies, and 1-dimensional torsional potentials at the `lvl_wbs` level and the energies at the `cc_lvl_d`. These are defined in the theory.dat file.

theory.dat::

    level wbs
        method = wb97xd
        basis = 6-31g*
        orb_res = RU
        program = gaussian09
    end level

    level ccdz
        method = ccsd(t)
        basis = cc-pvdz
        orb_res = RR
        program = molpro2015
    end level

Modify example for thermochem

Output
======

At the completion of ESDriver and kTPDriver, you will produce a MESS file and fit parameters.

MESS input file::

    MESS input

Note that fake wells have been added

CHEMKIN output::

    Rate params
