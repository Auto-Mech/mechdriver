""" Run OneDMin
"""

import automol.geom
import onedmin_io
import elstruct.writer
from autorun._run import from_parallel_input_strings


SCRIPT_NAME = 'run_onedmin.sh'
INPUT_NAME = 'input.dat'
OUTPUT_NAMES = ('output.dat', 'lj.out', 'min_geoms.dat', 'zero_ene')


# Specialized runners
def lennard_jones_params(sp_script_str, run_dir, nsamp, njobs,
                         tgt_geo, bath_geo, thy_info, charge, mult,
                         smin=2.0, smax=6.0, spin_method=1, ranseeds=None):
    """ Calculate an averaged Lennard-Jones epsilon and sigma parameter
        for the interaction potential between a target and bath molecule.

        :param sp_script_str: submission script for single-point calculation
        :type sp_script_str: str
        :param run_dir: directory where all OneDMin jobs are run
        :type run_dir: str
        :param nsamp: number of samples to run PER OneDMin job
        :type nsamp: int
        :param njobs: number of OneDMin instances to run in parallel
        :type njobs: int
        :param tgt_geo: geometry of the target molecule
        :type tgt_geo: automol geometry data structure
        :param bath_geo: geometry of the bath molecule
        :type bath_geo: automol geometry data structure
        :param thy_info: theory info object (prog, method, basis, orb_lbl)
        :type thy_info: tuple(str, str, str, str)
        :param charge: charge of the target-molecule complex
        :type charge: int
        :param mult: multiplicity of the target-molecule complex
        :type mult: int
        :param smin: minimum allowed intermolecular separation
        :type smin: float
        :param smax: maximum allowed intermolecular separation
        :type smax: float
        :param spin_method: parameter for the spin method
        :type spin_method: int
        :param ranseed: seed-integer for the orientational sampling
        :type ranseed: int
        :rtype: (float, float)
    """

    # maybe set the number of ranseeds to number of jobs?
    assert njobs == len(ranseeds)

    _, _, output_strs_lst = direct(
        sp_script_str, run_dir, nsamp, njobs,
        tgt_geo, bath_geo, thy_info, charge, mult,
        smin=smin, smax=smax, spin_method=spin_method, ranseeds=ranseeds)

    # Parse out the lj parameters and take the average
    sigmas, epsilons = (), ()
    for output_strs in output_strs_lst:
        sigmas, epsilons = onedmin_io.reader.lennard_jones(output_strs[1])
        sigmas += sigmas
        epsilons += epsilons

    return sigmas, epsilons


# General runners
def direct(sp_script_str, run_dir, nsamp, njobs,
           tgt_geo, bath_geo, thy_info, charge, mult,
           smin=3.779, smax=11.339, spin_method=1, ranseeds=None,
           script_name=SCRIPT_NAME,
           input_name=INPUT_NAME,
           output_names=OUTPUT_NAMES):
    """ Write input and run output.

        :param sp_script_str: submission script for single-point calculation
        :type sp_script_str: str
        :param run_dir: directory where all OneDMin jobs are run
        :type run_dir: str
        :param nsamp: number of samples to run PER OneDMin job
        :type nsamp: int
        :param njobs: number of OneDMin instances to run in parallel
        :type njobs: int
        :param tgt_geo: geometry of the target molecule
        :type tgt_geo: automol geometry data structure
        :param bath_geo: geometry of the bath molecule
        :type bath_geo: automol geometry data structure
        :param thy_info: theory info object (prog, method, basis, orb_lbl)
        :type thy_info: tuple(str, str, str, str)
        :param charge: charge of the target-molecule complex
        :type charge: int
        :param mult: multiplicity of the target-molecule complex
        :type mult: int
        :param smin: minimum allowed intermolecular separation
        :type smin: float
        :param smax: maximum allowed intermolecular separation
        :type smax: float
        :param spin_method: parameter for the spin method
        :type spin_method: int
        :param ranseed: seed-integer for the orientational sampling
        :type ranseed: int
        :rtype: (float, float)
    """

    # Write the main input files for all runs (breaks if ranseeds not given)
    input_str_lst = ()
    for ranseed in ranseeds:
        input_str_lst += (
            onedmin_io.writer.input_file(
                nsamp, smin, smax,
                ranseed=ranseed, spin_method=spin_method),
        )

    # Write the aux inputs; same for all runs
    tgt_str = automol.geom.string(tgt_geo)
    bath_str = automol.geom.string(bath_geo)

    elstruct_inp_str, onedmin_exe_name = _set_pot_info(thy_info, charge, mult)

    aux_dct = {
        'target.xyz': tgt_str,
        'bath.xyz': bath_str,
        'qc.mol': elstruct_inp_str,
        'qc.x': sp_script_str
    }

    # Write the script string for submission (for all runs)
    script_str = onedmin_io.writer.submission_script(
        njobs, run_dir, onedmin_exe_name)

    # Run the code
    output_str_lst = from_parallel_input_strings(
        script_str, run_dir, input_str_lst,
        aux_dct=aux_dct,
        script_name=script_name,
        input_name=input_name,
        output_names=output_names)

    return input_str_lst, elstruct_inp_str, output_str_lst


def _set_pot_info(thy_info, charge, mult):
    """ Figure out what the executables and elstruct should be based
        on the desired thy info.

        :param thy_info: theory info object (prog, method, basis, orb_lbl)
        :type thy_info: tuple(str, str, str, str)
        :param charge: charge of the target-molecule complex
        :type charge: int
        :param mult: multiplicity of the target-molecule complex
        :type mult: int
        :rtype: (str, str)
    """

    prog, method, basis, orb_lbl = thy_info

    if prog == 'exp6':
        elstruct_inp_str = None
        onedmin_exe_name = 'onedmin-exp6.x'
    else:
        elstruct_inp_str = elstruct.writer.energy(
            geo='GEOMETRY',
            charge=charge,
            mult=mult,
            method=method,
            basis=basis,
            prog=prog,
            mol_options=('nosym', 'noorient'),
            # memory=kwargs['memory'],
            memory=10.0,
            comment='SAMPLE GEOM',
            orb_type=orb_lbl
        )
        if 'gaussian' in prog:
            onedmin_exe_name = 'onedmin-dd-gaussian.x'
        else:
            onedmin_exe_name = 'onedmin-dd-molpro.x'

    return elstruct_inp_str, onedmin_exe_name
