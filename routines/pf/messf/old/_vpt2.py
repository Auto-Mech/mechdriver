"""
Functions to deal with xmat for projected freqs from QTC
"""

import os
import sys
import numpy as np

# from . import iotools as io
# from . import qctools as qc
import logging


def get_freqs(filename):
    """
    Pulls the frequencies out from EStokTP me output file 
    INPUT:
    filename - name of EStokTP output file (reac1_fr.me or reac1_unpfr.me)
    OUTPUT:
    freqs    - frequencies obtained from output file
    order    - in case the frequencies were reordered when sorting, keeps 
               track of which index of freqs corresponds to which normal mode
    """ 
    full = io.read_file(filename)
    full = full.strip('\n')
    full = full.split('[1/cm]')[1].split('Zero')[0] 
    full = full.split('ElectronicLevels')[0] 
    full = full.split('End')[0] 
    full = full.split()
    nfreqs = full[0]
    freqs = full[1:]
    if 'ts' in filename:
        imaginary = io.read_file(filename)
        imaginary=imaginary.split('ImaginaryFrequency[1/cm]')[1].split()[0]
        freqs.insert(0, '-' + imaginary)
    #[freq=float(freq) for freq in freqs]
    freqs = np.array(list(map(float, freqs)))
    a= freqs.argsort()
    freqs = np.sort(freqs)
    return freqs.tolist(), a.tolist()


def find_hinfreqs(proj, unproj, order):
    """
    Compares the frequencies from EStokTP projected and unprojected frequency
    output to determine which normal modes are hindered rotors
    INPUTS:
    proj   -  frequencies after projection
    unproj -  unprojected frequencies
    order  -  in case the frequencies were reordered when sorting, keeps track of 
              which index of unproj corresponds to which normal mode
    """
    diff = len(unproj) - len(proj)
    if diff > 0:
        for i in range(len(proj)):
            length = len(unproj)-1
            closeenough = 0.02
            for k in range(len(unproj)):
                if (abs(proj[i]-unproj[k]) < unproj[k] * closeenough):
                    #proj = np.delete(proj, 0)
                    unproj = np.delete(unproj, k)


def find_hinfreqs(proj, unproj, order):
    """
    Compares the frequencies from EStokTP projected and unprojected frequency
    output to determine which normal modes are hindered rotors
    INPUTS:
    proj   -  frequencies after projection
    unproj -  unprojected frequencies
    order  -  in case the frequencies were reordered when sorting, keeps track of 
              which index of unproj corresponds to which normal mode
    """
    diff = len(unproj) - len(proj)
    if diff > 0:
        for i in range(len(proj)):
            length = len(unproj)-1
            closeenough = 0.02
            for k in range(len(unproj)):
                if (abs(proj[i]-unproj[k]) < unproj[k] * closeenough):
                    #proj = np.delete(proj, 0)
                    unproj = np.delete(unproj, k)
                    order = np.delete(order, k)
                    break
        order = order[:diff]
    else:
        order = []
    modes = [mode+1 for mode in order]
    return modes


def remove_modes(xmat, modes):
    """
    Removes specified modes from anharmonic constant matrix
    INPUTS:
    xmat  - anharmonic constant matrix
    m?odes - the modes to delete from the matrix (with 1 being the first mode)
    OUTPUTS:
    xmat  - anharmonic constant matrix with columns and rows deleted for specified modes
    """
    modes.sort()#reverse=True)
    modeindex = [mode-1 for mode in modes]
    for index in modeindex[::-1]:
        xmat = np.delete(xmat, index, 0)
        xmat = np.delete(xmat, index, 1)
    return xmat


def remove_vibrots(vibrot, modes):
    """
    Removes specified modes from anharmonic constant matrix
    INPUTS:
    xmat  - anharmonic constant matrix
    modes - the modes to delete from the matrix (with 1 being the first mode)
    OUTPUTS:
    xmat  - anharmonic constant matrix with columns and rows deleted for specified modes
    """
    modes.sort()#reverse=True)
    vibrot = vibrot.splitlines()
    modeindex = [mode-1 for mode in modes]
    vibrots = []
    for index in range(len(vibrot)):
        if index not in modeindex:
            vibrots.append(vibrot[index])
    return '\n'.join(vibrots)


def anharm_freq(freqs, xmat):
    """
    Uses anharmonic frequency matrix and harmonic frequencies to compute VPT2 anharmonic frequencies
    INPUT:
    freqs   - harmonic frequencies
    xmat    - anharmonic constant matrix
    OUTPUT:
    anharms - VPT2 anharmonic frequencies
    """
    anharms = np.zeros(len(freqs))
    for i, freq in enumerate(freqs):
        anharms[i]  = freq
        anharms[i] += 2. * xmat[i][i]
        tmp = 0
        for j in range(len(freqs)):
            if j != i:
                tmp += xmat[i][j]
        if tmp > 0:
            logging.warning('Positive anharmonic correction on Mode {:d}'.format(i+1))
        if tmp < 400:
            logging.warning('Large anharmonic correction of {:f} on Mode {:d}'.format(tmp, i+1))
        anharms[i] += 1./2 * tmp

    return anharms


def anharm_zpve():
    """ Calculate the anharmonic zpve using the x mat
    """
    return anharm_zpve
# def main(args, vibrots = None):
#
#     if isinstance(args, dict):
#         if 'pfreqs' in args:
#             proj = np.array(args['pfreqs']).astype(np.float)
#             unproj = np.array(args['freqs']).astype(np.float)
#             proj = np.sort(proj)
#             unproj = np.sort(unproj)
#             a = np.arange(len(proj)+1)
#             b = np.arange(len(unproj)+1)
#             #a = nargs['pfreqs']).argsort()[::-1]
#             #b = args['freqs'].argsort()[::-1]
#         else:
#             proj, a   = get_freqs(args['freqfile'])
#             unproj, b = get_freqs(args['unprojfreq'])
#
#
#         xmat      = remove_modes(xmat, modes)
#         anfreq = anharm_freq(proj, xmat)
#         else:
#             xmat = []
#             anfreq = proj
#
#         if args.computeanharm.lower() == 'true':
#             xmat = gauss_xmat(anharmlog, natoms)
#             proj, b   = get_freqs(eskproj)
#             unproj, a = get_freqs(eskunproj)
#             modes     = find_hinfreqs(proj, unproj, a)
#             xmat      = remove_modes(xmat, modes)
#             proj, b   = get_freqs(eskproj)
#             anfreq = anharm_freq(proj, xmat)
