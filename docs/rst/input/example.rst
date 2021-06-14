
example
=======

Input
-----

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

