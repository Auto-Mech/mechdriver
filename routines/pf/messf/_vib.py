"""
  Handle vibrational data info
"""

import os
import projrot_io
import autofile

# New libs
from lib.phydat import phycon
from lib.runner import script


def projrot_freqs_1(tors_geo, hess,
                    proj_rotors_str,
                    save_path, pot=False, saddle=False):
    """ Get frequencies from one version of ProjRot
    """
    coord_proj = 'cartesian'
    grad = ''
    # Write the string for the ProjRot input
    projrot_inp_str = projrot_io.writer.rpht_input(
        tors_geo, grad, hess, rotors_str=proj_rotors_str,
        coord_proj=coord_proj)

    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(save_path)
    bld_save_fs.leaf.create(bld_locs)
    path = bld_save_fs.leaf.path(bld_locs)
    print('Build Path for Partition Functions in species block')
    print(path)
    proj_file_path = os.path.join(path, 'RPHt_input_data.dat')
    with open(proj_file_path, 'w') as proj_file:
        proj_file.write(projrot_inp_str)

    script.run_script(script.PROJROT, path)

    freqs = []
    zpe_har_no_tors = 0.
    # har_zpe = 0.
    if pot:
        rthrproj_freqs, _ = projrot_io.reader.rpht_output(
            path+'/hrproj_freq.dat')
        freqs = rthrproj_freqs
        zpe_har_no_tors = sum(freqs)*phycon.WAVEN2KCAL/2.
    rtproj_freqs, imag_freq = projrot_io.reader.rpht_output(
        path+'/RTproj_freq.dat')
    # har_zpe = sum(rtproj_freqs)*phycon.WAVEN2KCAL/2.
    if not freqs:
        freqs = rtproj_freqs
    if saddle:
        if imag_freq:
            imag_freq = imag_freq[0]
        else:
            imag_freq = freqs[-1]
            freqs = freqs[:-1]

    return freqs, imag_freq, zpe_har_no_tors


def projrot_freqs_2(save_path, pot=False, saddle=False):
    """ Get ProjRot frequencies via ProjRot 2
    """

    bld_locs = ['PROJROT', 0]
    bld_save_fs = autofile.fs.build(save_path)
    bld_save_fs.leaf.create(bld_locs)
    path = bld_save_fs.leaf.path(bld_locs)

    projrot_script_str2 = (
        "#!/usr/bin/env bash\n"
        "RPHt.exe >& /dev/null")
    script.run_script(projrot_script_str2, path)

    zpe_har_no_tors_2 = 0.0
    freqs_2 = []
    if pot:
        rthrproj_freqs_2, _ = projrot_io.reader.rpht_output(
            path+'/hrproj_freq.dat')
        freqs_2 = rthrproj_freqs_2
        zpe_har_no_tors_2 = sum(freqs_2)*phycon.WAVEN2KCAL/2.
    rtproj_freqs, imag_freq_2 = projrot_io.reader.rpht_output(
        path+'/RTproj_freq.dat')
    har_zpe = sum(rtproj_freqs)*phycon.WAVEN2KCAL/2.
    if not freqs_2:
        freqs_2 = rtproj_freqs
    if saddle:
        if imag_freq_2:
            imag_freq_2 = imag_freq_2[0]
        else:
            imag_freq_2 = freqs_2[-1]
            freqs_2 = freqs_2[:-1]

    return freqs_2, imag_freq_2, har_zpe, zpe_har_no_tors_2


def determine_freqs_zpe(freqs1, freqs2, imag_freq1, imag_freq2,
                        zpe_harm_no_tors, zpe_harm_no_tors_2,
                        harm_zpe, tors_zpe):
    """ get the freqs ftom two methods
    """
    harm_tors_zpe = harm_zpe - zpe_harm_no_tors
    harm_tors_zpe_2 = harm_zpe - zpe_harm_no_tors_2
    del_tors_zpe = harm_tors_zpe - tors_zpe
    del_tors_zpe_2 = harm_tors_zpe_2 - tors_zpe
    if del_tors_zpe <= del_tors_zpe_2:
        zpe = zpe_harm_no_tors + tors_zpe
        freqs = freqs1
        imag_freq = imag_freq1
    else:
        zpe = zpe_harm_no_tors_2 + tors_zpe
        freqs = freqs2
        imag_freq = imag_freq2
    if abs(del_tors_zpe) > 0.2 and abs(del_tors_zpe_2) > 0.2:
        print('Warning: There is a difference of ',
              '{0:.2f} and {1:.2f}'.format(del_tors_zpe, del_tors_zpe_2),
              'kcal/mol between harmonic and hindered torsional ZPVEs')

    return freqs, imag_freq, zpe
