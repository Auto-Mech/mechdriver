
Complex Example
---------------

Workflow
~~~~~~~~

The following are a set of inputs for a minimal example treatment of the C3H8+OH reaction system. This system consists of an initial hydrogen abstraction via OH, followed by decomposition of the product C3H7 radical. While not the simplest example, it does highlight several powerful features of the code.

This example was written to be used with the Psi4 program as an open source.

The following input will work on a clean system and does not assume any existing run-save filesystem. All required calculations should be performed. Although the same output should be achieved if run a second time.

The only modification necessary is in the run.dat file where the run_prefix and save_prefix paths. 

A simple theoretical treatment has been employed. Quantitative accuracy should not be expected.
    (1) electronic structure: low-level DFT methods, small basis-set MP2 energies
    (2) thermochemistry: RRHO(1DHR)
    (3) kinetics: RRHO(1DHR), with fixed transition state theory, Eckart tunneling; energy transfer params estimated with internal scheme

The steps of the workflow as follows:
    (1) ESDriver: Generate geometries, frequencies, and energies PES specified in pes block (C2H6, C2H5, H, H2, TS)
    (2) ESDriver: Generate geometries, frequencies, and energies PES specified in spc block (CH4)
    (3) ThermoDriver: Build and run MESSPF for partition functions, then generate NASA polynomials for all species 
    (4) kTPDriver: Build and run MESSRATEs for rate constants, then fit them


Input
~~~~~

To set up the chemical reactions and species for the input mechanism, we set

mechanism.dat file::

    REACTIONS
        C3H8+OH=C3H7(1)+H2O  1.0  0.0  0.0
        C3H8+OH=C3H7(2)+H2O  1.0  0.0  0.0
        C3H7(1)=C3H7(2)      1.0  0.0  0.0
        C3H7(1)=C3H6+H       1.0  0.0  0.0
        C3H7(2)=C3H6+H       1.0  0.0  0.0
    END

species.csv file::

    name,smiles,mult,charge
    C3H8,'CCC',1,0
    C3H7(1),'CC[CH2]',2,0
    C3H7(2),'C[CH]C',2,0
    C3H6(1),'CC=C',1,0
    OH,'[OH]',2,0
    H,'[H]',2,0
    H2O,'[HH]',1,0
    CH4,'CC',1,0

run.dat file::

    input
        run_prefix = /fake/path/to/run
        save_prefix = /fake/path/to/save
    end input

    pes
        1: 1
    end pes

    spc
        8
    end spc

    els
        spc init_geom     runlvl=wbsgs   inplvl=wbsgs
        spc conf_samp     runlvl=wbsgs   inplvl=wbsgs
        ts  find_ts       runlvl=wbsgs   inplvl=wbsgs
        all hr_scan       runlvl=wbsgs   inplvl=wbsgs  tors_model=1dhr
        all conf_energy   runlvl=mp2dz   inplvl=wbsgs
        all conf_hess     runlvl=wbsgs   inplvl=wbsgs
    end els

    thermo
        write_mess      kin_model=global  spc_model=global
        run_mess        kin_model=global  spc_model=global
        run_fits        kin_model=global  spc_model=global
    end thermo

    ktp
        write_mess      kin_model=global  spc_model=global
        run_mess
        run_fits        kin_model=global  spc_model=global
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

    level wbsgs
        method = b3lyp
        basis = 6-31g*
        orb_res = RU
        program = psi4
    end level

    level mp2dz
        method = mp2
        basis = cc-pvdz
        orb_res = RR
        program = psi4
    end level

Modify example for thermochem

Output
~~~~~~

At the completion of ESDriver and kTPDriver, you will produce a MESS file and fit parameters.

MESS input file::

    MESS input STR

Note that fake wells have been added

CHEMKIN output::

    Rate params

