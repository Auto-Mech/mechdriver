""" moldr drivers
"""
import autodir
import elstruct
import moldr.run
import moldr.save


def conformers(nsamp, script_str, run_prefix, save_prefix,
               # elstruct robust run arguments
               prog, method, basis, geo, mult, charge,
               errors=(), options_mat=(),
               **kwargs):
    """ search for unique conformers
    """
    nsamp_tot = nsamp

    if autodir.conf.has_base_information_file(save_prefix):
        inf_obj = autodir.conf.read_base_information_file(save_prefix)
        nsamp = max(nsamp_tot - inf_obj.nsamp, 0)
        nsamp_tot = max(nsamp_tot, inf_obj.nsamp)
        print("Found previous saved run. Adjusting `nsamp`.")
        print("nsamp done={:d}".format(inf_obj.nsamp))
        print("nsamp left={:d}".format(nsamp))

    moldr.run.torsion_sampling(
        nsamp=nsamp,
        script_str=script_str,
        prefix=run_prefix,
        input_writer=elstruct.writer.optimization,
        prog=prog,
        method=method,
        basis=basis,
        geo=geo,
        mult=mult,
        charge=charge,
        errors=errors,
        options_mat=options_mat,
        **kwargs,
    )

    moldr.save.conformers(
        run_prefix=run_prefix,
        save_prefix=save_prefix,
    )
