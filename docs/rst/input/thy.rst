
theory.dat
----------

Details various electronic structure methods.

Theory levels are given as a series of blocks like so:

General format::

    level <level_name>
            <keyword1 = value>
            <keyword2 = value>
              …
    end level

Above, <level_name> serves as a descriptor of the method that is used in other parts of the code. This string can be whatever the user wishes. It must contain no whitespace.

Example::

    level mp2gs
        method = mp2
        basis = 6-31g*
        orb_res = RU
        program = psi4
    end level

The following example shows the minimum required for a `fully defined` block. Fully
defines the method and basis. Set how the references are set based on closed/open shell
treatments. Sets program.

However additional keywords may be given to specify either more options or job-time params. If jobtime params (like nprocs, memorY) are not set, we set these values according
to the program.

More Example of two blocks with same method but different runtime options (core-valence)::

    level mp2gs
        method = mp2
        basis = 6-31g*
        orb_res = RU
        program = psi4
        memory = 10
        nprocs = 8
    end level

Here we have modifed the method to have runtime options

More Example of two blocks with same method but different runtime options (core-valence)::

    level fc_ccdz_psi4
        method = ccsd(t)
        basis = cc-pcvdz
        orb_res = RU
        program = psi4
    end level

    level ae_ccdz_psi4
        method = ccsd(t)
        basis = cc-pcvdz
        orb_res = RU
        program = psi4
        job_options = [all_elec]
    end level

See that these two blocks are largely similar, only differing in the job options.

The range of capabilities forallowed methods is enforced by the elstruct package.

Any number of these level blocks can be given in theory.dat. The only blocks used by the code are those specified in run.dat and models.dat. (Even though all of the blocks will be read)

IMPORTANT: methods currently treated as consistent across all programs.

