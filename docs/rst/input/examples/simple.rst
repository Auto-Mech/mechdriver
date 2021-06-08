
Simple Example
--------------

Workflow
~~~~~~~~

Following for a very simple example of the C2H6+H reaction system. 

A simple theoretical treatment has been employed. Quantitative accuracy should not be expected.
    (1) electronic structure: low-level DFT methods, small basis-set MP2 energies
    (2) thermochemistry: RRHO
    (3) kinetics: RRHO with fixed transition state theory; energy transfer params estimated with internal scheme

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
        C2H6+H=C2H5+H2  1.0  0.0  0.0
    END

species.csv file::

    name,smiles,mult,charge
    C2H6,'CC',1,0
    C2H5,'C[CH2]',2,0
    H,'[H]',2,0
    H2,'[HH]',1,0
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
        5
    end spc

    els
        spc init_geom     runlvl=wbsgs   inplvl=wbsgs
        ts  find_ts       runlvl=wbsgs   inplvl=wbsgs
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

