""" temporary script
"""
import elstruct


def first_successful_output(script_str, elstruct_input_writer,
                            prog, method, basis, geom, mult, charge,
                            try_kwargs_lst, **kwargs):
    """ run through several options and return the first succesful output
    """
    for try_kwargs in try_kwargs_lst:
        assert not set(try_kwargs.keys()) & set(kwargs.keys())
        inp_str = elstruct_input_writer(
            prog=prog, method=method, basis=basis, geom=geom, mult=mult,
            charge=charge, **try_kwargs, **kwargs)
        out_str, _ = elstruct.run(script_str, inp_str)
        if elstruct.reader.ran_successfully(prog, out_str):
            break
    return out_str


if __name__ == '__main__':
    import autofile

    # input arguments
    SCRIPT_STR = "#!/usr/bin/env bash\npsi4 >> stdout.log &> stderr.log"
    PROG = 'psi4'
    METHOD = 'rhf-mp2'
    BASIS = 'sto-3g'
    GEOM = (('O', (None, None, None), (None, None, None)),
            ('H', (0, None, None), (1.8, None, None)),
            ('H', (0, 1, None), (1.8, 1.75, None)))
    MULT = 1
    CHARGE = 0
    ZMAT_VAR_DCT = {(1, 0): 'r1', (2, 0): 'r1'}

    TRY_KWARGS_LST = [
        {'scf_options': 'set guess core'},
        {'scf_options': 'set guess gwh'},
        {'scf_options': 'set guess sad'},
    ]

    # runner
    OUT_STR = first_successful_output(
        SCRIPT_STR, elstruct.writer.optimization,
        PROG, METHOD, BASIS, GEOM, MULT, CHARGE,
        TRY_KWARGS_LST, zmat_var_dct=ZMAT_VAR_DCT)

    # call readers to get information from output
    ENE = elstruct.reader.energy(PROG, METHOD, OUT_STR)
    GEO = elstruct.reader.optimized_geometry(PROG, OUT_STR)
    RET = elstruct.reader.optimized_zmatrix(PROG, OUT_STR)

    # write things to files
    autofile.write.energy('h2o', ENE)
    autofile.write.geometry('h2o', GEO)
    if RET is not None:
        ZMA, VAR_DCT = RET
        autofile.write.zmatrix('h2o', ZMA, var_dct=VAR_DCT)
