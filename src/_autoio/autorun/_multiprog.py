""" Combined process for getting frequencies
"""

import os
import automol.form
from phydat import phycon
import mess_io.reader
import thermp_io.reader
import projrot_io.writer
from autorun.projrot import frequencies
from autorun.thermp import direct as thermp_direct
from autorun.pac99 import nasa_polynomial
from autorun.mess import torsions as mess_torsions


# PROJROT + MESS Runners
def projected_frequencies(mess_script_str, projrot_script_str, run_dir,
                          mess_hr_str, projrot_hr_str,
                          mess_geo, projrot_geo, hess,
                          dist_cutoff_dct1=None, dist_cutoff_dct2=None,
                          saddle=False):
    """ Compute the harmonic vibrational frequencies and ZPVE after projecting
        out the hindered rotors. Tests different cutoffs for defininig rotors
        and sees which reproduces the harmonic frequency result; takes that one
    """

    # Calculate the torsional frequencies using MESS
    tors_freqs, _ = mess_torsions(
        mess_script_str, run_dir, mess_geo, mess_hr_str)

    # Calculate the projected vibrational frequencies using ProjRot
    if dist_cutoff_dct1 is None:
        dist_cutoff_dct1 = {('H', 'O'): 2.26767, ('H', 'C'): 2.26767}
    rotor_dist1_str = projrot_io.writer.projection_distance_aux(
        dist_cutoff_dct=dist_cutoff_dct1)
    aux_dct1 = {'dist_rotpr.dat': rotor_dist1_str}

    print('running projrot the first time:')
    run_dir_1 = os.path.join(run_dir, '1')
    rt_freqs1, rth_freqs1, rt_imag1, _ = frequencies(
        projrot_script_str, run_dir_1, [projrot_geo], [[]], [hess],
        rotors_str=projrot_hr_str, aux_dct=aux_dct1)

    if dist_cutoff_dct2 is None:
        dist_cutoff_dct2 = {('H', 'O'): 2.83459, ('H', 'C'): 2.83459,
                            ('C', 'O'): 3.7807}
    rotor_dist2_str = projrot_io.writer.projection_distance_aux(
        dist_cutoff_dct=dist_cutoff_dct2)
    aux_dct2 = {'dist_rotpr.dat': rotor_dist2_str}

    print('running projrot the second time:')
    run_dir_2 = os.path.join(run_dir, '2')
    _, rth_freqs2, rt_imag2, _ = frequencies(
        projrot_script_str, run_dir_2, [projrot_geo], [[]], [hess],
        rotors_str=projrot_hr_str, aux_dct=aux_dct2)

    if not rth_freqs2:
        dist_cutoff_dct3 = {('H', 'O'): 3.401506, ('H', 'C'): 3.779451,
                                ('C', 'O'): 4.53534}
        rotor_dist3_str = projrot_io.writer.projection_distance_aux(
            dist_cutoff_dct=dist_cutoff_dct3)
        aux_dct3 = {'dist_rotpr.dat': rotor_dist3_str}

        print('running projrot a third time:')
        run_dir_3 = os.path.join(run_dir, '3')
        _, rth_freqs2, rt_imag2, _ = frequencies(
            projrot_script_str, run_dir_3, [projrot_geo], [[]], [hess],
            rotors_str=projrot_hr_str, aux_dct=aux_dct3)
        if rth_freqs2:
            print(
                'it did work, make sure these frequencies look alright:',
                rth_freqs2)
    # Calculate ZPVEs from all harmonic freqs and torsional freqs
    tors_zpe = (sum(tors_freqs) / 2.0) * phycon.WAVEN2EH
    harm_zpe = (sum(rt_freqs1) / 2.0) * phycon.WAVEN2EH

    # Calculate harmonic ZPVE from freqs where torsions have been projected out
    # Value from both projrot versions, which use different projection schemes
    harm_zpe_notors_1 = (sum(rth_freqs1) / 2.0) * phycon.WAVEN2EH
    harm_zpe_notors_2 = (sum(rth_freqs2) / 2.0) * phycon.WAVEN2EH

    # Calcuate the difference in the harmonic ZPVE from projecting out torsions
    harm_tors_zpe = harm_zpe - harm_zpe_notors_1
    harm_tors_zpe_2 = harm_zpe - harm_zpe_notors_2

    # Check to see which of the above ZPVEs match more closely with tors ZPVE
    # calculated directly by treating the torsions in MESS
    diff_tors_zpe = harm_tors_zpe - tors_zpe
    diff_tors_zpe_2 = harm_tors_zpe_2 - tors_zpe
    if diff_tors_zpe <= diff_tors_zpe_2:
        proj_freqs = rth_freqs1
        proj_imag = rt_imag1
        proj_zpe = harm_zpe_notors_1
    else:
        proj_freqs = rth_freqs2
        proj_imag = rt_imag2
        proj_zpe = harm_zpe_notors_2

    # Check imaginary frequencies and set freqs
    if saddle:
        if len(proj_imag) > 1:
            print(
                'There is more than one imaginary frequency')
        proj_imag = max(proj_imag)
    else:
        for _ in proj_imag:
            rt_freqs1 += (0.00001,)
        proj_imag = None

    # Check if there are significant differences caused by the rotor projection
    diff_tors_zpe *= phycon.EH2KCAL
    diff_tors_zpe_2 *= phycon.EH2KCAL
    if abs(diff_tors_zpe) > 0.2 and abs(diff_tors_zpe_2) > 0.2:
        print(
            'The first and second runs yield differences of ',
            f'{diff_tors_zpe:.2f} and {diff_tors_zpe_2:.2f}',
            'kcal/mol between harmonic and hindered torsional ZPVEs')

    return proj_freqs, proj_imag, proj_zpe, rt_freqs1, tors_freqs


# ThermP + PAC99 Runners
def thermo(thermp_script_str, pac99_script_str, run_dir,
           pf_str, name, formula, hform0,
           enthalpyt=0.0, breakt=1000.0, convert=False):
    """ Runs ThermP and PAC99 together and parse out the results

        return hform298 k in Eh
    """

    # Obtain the temperatures from the pf.dat file used to fit the thermo
    temps, _, _, _ = mess_io.reader.pfs.partition_function(pf_str)
    temps = temps[:-1]

    # Run ThermP and read the Heat-of-Formation at 298K
    formula_str = automol.form.string(formula)
    _, thermp_output_strs = thermp_direct(
        thermp_script_str, run_dir,
        pf_str, formula_str, hform0, temps,
        enthalpyt=enthalpyt, breakt=breakt)
    hform298 = thermp_io.reader.hf298k(thermp_output_strs[0])

    # Run PACC99 and obtain the polynomials
    nasa_poly = nasa_polynomial(
        pac99_script_str, run_dir, thermp_output_strs[1], name, formula,
        convert=convert)

    return hform298, nasa_poly
