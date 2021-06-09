
models.dat
----------

Easily the most complicated input file as here we define the mixtures of all the electronic structure methods as well as assumptions and treatments needed to build a full model of a species for thermo calcs and rate calculations.


Only required for thermo, ktp, and proc driver.

Two types of sections can be defined in the file::

    kin <model_name>
        <model information>
    end <model_name>
    
    spc <model_name>
        <model information>
    end <model_name>

The kin block is used to define conditions and base assumptions used to calculate and fit kinetics and thermochemistry for entire PES or species lists.

The spc block is used to define the assumptions used to calculate the individual partition functions for species and transition states comprising PES or species list.

Each block consists of several parenthetical keyword blocks which may consist lists of numbers or keyword-value pairs:

kin example::

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



spc example::

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

Each subsection defines various keywords for evaluating partitions. For most sections (ene, rot, vib) there are `mod` and `lvl` keywords.

Arranged to build each section using any general combination of methods and assumptions to calculate the pfs.

    mod: central model/assumption for evaluating part'n fxn.
    lvl: electronic structure method(s) used for data. must correspond to data that exists in save filesystem. may be fully defined by series/combinations of methods
    geolvl: geometry method
    enelvl: energy for single points
    lvln: sequence of lvls for composite methods

For ts there is additional considerations for transition state theory methods.

As with many other files any number of such named blocks may be used.

