""" Handle arithmetic operations for partition functions and
    thermodynamic quantities derived from them, including
    enthalpy, entropy, gibbs, etc
"""

import numpy

from phydat import phycon
import automol.geom
# import autorun
# from autorun.thermp import direct as thermp_direct
# import thermp_io
# import mess_io


def combine(pfs, coeffs, operators):
    """ Combine sets of partition function data together via some set
        of input coefficients and arithmetic operations.

        _pf = [temps, logq, dq_dt, d2q_dt2]
    """

    assert len(pfs) == len(coeffs)

    for midx, _pf in enumerate(pfs):
        if midx == 0:
            coeff = coeffs[midx]
            final_pf = _pf
        else:
            pf2 = _pf
            coeff = coeffs[midx]
            operator = operators[midx-1]
            if coeff < 0:
                coeff = abs(coeff)
                if operator == 'multiply':
                    operator = 'divide'
            print(pf2)
            final_pf = _combine_pfs(final_pf, pf2, coeff, operator)

    return final_pf


def _combine_pfs(pfa, pfb, coeff, operator):
    """ Combine two sets of partition function data together under
        multiplication or division where the combined data is
        weighted by some coefficient.

        pfn = ((temp1, logq1, dqdt1, d2gdt21), (temp2, logq2, dqdt2, d2gdt22),)

        :param pfa: partition function set A
        :param pfb: partition function set B
        :param coeff: weighting coefficient
        :param operator: multiplication/division operator
    """

    tempsa, logqa, dq_dta, d2q_dt2a = pfa
    _, logqb, dq_dtb, d2q_dt2b = pfb
    if operator == 'multiply':
        logq = [a+b+numpy.log(coeff) for a, b in zip(logqa, logqb)]
        dq_dt = [a+b+numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
        d2q_dt2 = [a+b+numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]
    elif operator == 'divide':
        logq = [a-b-numpy.log(coeff) for a, b in zip(logqa, logqb)]
        dq_dt = [a-b-numpy.log(coeff) for a, b in zip(dq_dta, dq_dtb)]
        d2q_dt2 = [a-b-numpy.log(coeff) for a, b in zip(d2q_dt2a, d2q_dt2b)]

    return tempsa, tuple(logq), tuple(dq_dt), tuple(d2q_dt2)


def q_translational(mass, temp):
    """ Caulculate the translational partition function
    """
    return ((2 * numpy.pi * (mass * phycon.AMU2KG) * phycon.KB * temp)**(3/2)
            * (phycon.H * 100)**(-3))


def q_rotational(i_a, i_b, i_c, sigma, temp, linear=False):
    """ Caulculate the rotational partition function
    """

    def _aval(i_x):
        """ Convert moment of inertia to rot constant
        """
        return phycon.H / (8*numpy.pi**2 * i_x)

    def _q_rotational_linear(i_a, i_b, i_c, sigma, temp):
        """ Calculate linear term and convert from amu*bohr^2 in kg*m^2
        """
        i_a = max(i_a, i_b, i_c)
        i_a = i_a * phycon.AMU2KG * phycon.BOHR2CM**2 / 100**2
        return (phycon.KB * temp) / (phycon.H * _aval(i_a) * sigma)

    def _q_rotational_nonlinear(i_a, i_b, i_c, sigma, temp):
        """ Calculate nonlinear term and convert from amu*bohr^2 in kg*m^2
        """
        i_a = i_a * phycon.AMU2KG * phycon.BOHR2CM**2 / 100**2
        i_b = i_b * phycon.AMU2KG * phycon.BOHR2CM**2 / 100**2
        i_c = i_c * phycon.AMU2KG * phycon.BOHR2CM**2 / 100**2
        return ((phycon.KB * temp / phycon.H)**(3/2) * numpy.pi**(1/2)
                * (_aval(i_a) * _aval(i_b) * _aval(i_c))**(-1/2) / sigma)

    if linear:
        q_rot = _q_rotational_linear(i_a, i_b, i_c, sigma, temp)
    else:
        q_rot = _q_rotational_nonlinear(i_a, i_b, i_c, sigma, temp)
    return q_rot


def q_vibrational(freqs, temp):
    """ Caulculate the vibrational partition function
    """
    q_vib = 1.
    for freq in freqs:
        # nu_i = freq * phycon.WAVEN2EH * phycon.EH2KJ * 1000.  # in J/mol
        nu_i = freq * phycon.SOLMS * 100.   # in 1/s
        # numerator = numpy.exp(-(phycon.H * nu_i) / (2 * phycon.KB * temp))
        numerator = 1.
        denominator = 1 - numpy.exp(-(phycon.H * nu_i) / (phycon.KB * temp))
        q_vib *= numerator / denominator
    return q_vib


def rrho_partition_function(geo, freqs, temp_range=None, nlog=0):
    """ Calculate the total rrho partition function
    """
    q_total = {}
    if temp_range is None:
        temp_range = list(range(300, 3000, 100))
    mass = automol.geom.total_mass(geo)
    moms = automol.geom.moments_of_inertia(geo)
    linear = automol.geom.is_linear(geo)
    ext_symm = automol.geom.external_symmetry_factor(geo)
    int_symm, _ = automol.symm.oxygenated_hydrocarbon_symm_num(geo)
    sigma = int_symm * ext_symm
    for temp in temp_range:
        q_elec = 1.
        q_trans = q_translational(mass, temp)
        q_rot = q_rotational(*moms, sigma, temp, linear=linear)
        q_vib = q_vibrational(freqs, temp)
        q_all = q_elec * q_trans * q_rot * q_vib
        if nlog == 0:
            q_total[round(temp, 4)] = q_all
        elif nlog == 1:
            q_total[round(temp, 4)] = numpy.log(q_all)
        elif nlog == 2:
            q_total[round(numpy.log(temp), 4)] = numpy.log(q_all)
    return q_total


def pf_polys(pf_temp_dct, order=3):
    """ Fit partition function to polynomial and take the derivatives?
    """
    pf_array, temp_array = list(pf_temp_dct.values()), list(pf_temp_dct.keys())
    poly = numpy.polyfit(temp_array, pf_array, order)
    func = numpy.poly1d(poly)
    der = func.deriv()
    der2 = der.deriv()
    return func, der, der2


def heat_capacity_from_pf(lnq, dlnqdt, d2lnqdt2, temp):
    """ Calculate the heat capacity from the partition functon. [units?]
    """
    _ = lnq  # to get by pylint, remove lnq from function args???
    return (
        phycon.NAVO * phycon.KB
        * (temp**2 * d2lnqdt2(temp) + 2 * dlnqdt(numpy.log(temp)) + 1)
        ) * phycon.J2CAL


def enthalpy_from_pf(pf_fun, dqdt, temp):
    """ Calculate the enthalpy from the partition functon. [units?]
    """
    return phycon.RC_KCAL * temp * (temp / pf_fun(temp) * dqdt(temp) + 1)


def entropy_from_pf(pf_fun, dqdt, temp):
    """ WORKS with mystical 5.4 number
        Calculate the entropy from the partition functon. [units?]
    """
    entropy = phycon.NAVO * phycon.KB * (
        temp/pf_fun(temp)*dqdt(temp)
        + numpy.log(pf_fun(temp))
        - numpy.log(phycon.NAVO) + 5.4
        # - numpy.log(phycon.NAVO) + 1
        + numpy.log(temp))
    return entropy * phycon.J2CAL


def gibbs_energy_from_pf(pf_fun, temp):
    """ Calculate the internal energy from the partition functon. [units?]
    """
    energy = - phycon.NAVO * phycon.KB * temp * numpy.log(pf_fun(temp))
    return energy * phycon.J2CAL / 1000.


def rel_gibbs_energy_from_pf(pf_fun, temp, zero_ene):
    """ Calcluate the relative Gibbs energy from pf
    """
    rel_pf = (
        pf_fun(temp) * numpy.exp(-zero_ene / (phycon.KB * phycon.NAVO * temp))
    )
    energy = - phycon.NAVO * phycon.KB * temp * numpy.log(rel_pf)
    return energy * phycon.J2CAL / 1000.


def rrho_del_enthalpy(geo, freqs, temp=298.15):
    """ Get enthalpy from RRHO?
    """
    temp_range = numpy.arange(temp, temp+20, .05)
    q_total = rrho_partition_function(geo, freqs, temp_range, nlog=0)
    pf_fun, dqdt, _ = pf_polys(q_total)
    enthalpy = enthalpy_from_pf(pf_fun, dqdt, temp)
    return enthalpy


def rrho_entropy(geo, freqs, temp=298.15):
    """ Entropy from RRHO?
    """
    temp_range = numpy.arange(temp, temp+20, .05)
    q_total = rrho_partition_function(geo, freqs, temp_range, nlog=0)
    pf_fun, dqdt, _ = pf_polys(q_total)
    entropy = entropy_from_pf(pf_fun, dqdt, temp)
    return entropy


def rrho_heat_capacity(geo, freqs, temp=298.15):
    """ Heat Capacity from RRHO?
    """
    temp_range = numpy.arange(temp, temp+20, .05)
    heat_cap = None
    if temp > 20:
        lnq_total = rrho_partition_function(geo, freqs, temp_range, nlog=1)
        lnq, _, d2lnqdt2 = pf_polys(lnq_total)
        lnq_lnt = rrho_partition_function(geo, freqs, temp_range, nlog=2)
        _, dlnqdlnt, _ = pf_polys(lnq_lnt)
        heat_cap = heat_capacity_from_pf(lnq, dlnqdlnt, d2lnqdt2, temp)
    return heat_cap


def rrho_gibbs(geo, freqs, temp=298.15):
    """ RRHO Gibbs function
    """
    temp_range = numpy.arange(temp, temp+20, .05)
    q_total = rrho_partition_function(geo, freqs, temp_range, nlog=0)
    pf_fun, _, _ = pf_polys(q_total)
    gibbs = gibbs_energy_from_pf(pf_fun, temp)
    return gibbs


def rrho_gibbs_factor(geo, freqs, zero_ene, temp):
    """ RRHO Gibbs factor
    """
    zero_ene = zero_ene * 1000. / phycon.J2CAL
    temp_range = numpy.arange(temp, temp+20, .05)
    q_total = rrho_partition_function(geo, freqs, temp_range, nlog=0)
    pf_fun, _, _ = pf_polys(q_total)
    return rel_gibbs_energy_from_pf(pf_fun, temp, zero_ene)


def rrho_properties(geo, freqs, temps=None):
    """ RRHO Props
    """
    if temps is None:
        temps = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1500]
    for temp in temps:
        temp_range = numpy.arange(temp, temp+20, .05)
        heat_cap = 0
        if temp > 20:
            lnq_total = rrho_partition_function(geo, freqs, temp_range, nlog=1)
            lnq, _, d2lnqdt2 = pf_polys(lnq_total)
            # print('Q:', temp, lnq(temp), dlnqdt(temp), d2lnqdt2(temp))
            lnq_lnt = rrho_partition_function(geo, freqs, temp_range, nlog=2)
            _, dlnqdlnt, _ = pf_polys(lnq_lnt)
            heat_cap = heat_capacity_from_pf(lnq, dlnqdlnt, d2lnqdt2, temp)
        q_total = rrho_partition_function(geo, freqs, temp_range, nlog=0)
        pf_fun, dqdt, _ = pf_polys(q_total)
        entropy = entropy_from_pf(pf_fun, dqdt, temp)
        enthalpy = enthalpy_from_pf(pf_fun, dqdt, temp)
        # gibbs = enthalpy - entropy * temp / 1000.
        print('Prop:', temp, heat_cap, entropy, enthalpy)
    return enthalpy, entropy, heat_cap


# def gibbs_from_property_dct(csh_t_dct, hform0k):
#     """
#     """
#     gibbs_t_dct = {}
#     for temp in csh_t_dct:
#         _, entropy, enthalpy = csh_t_dct[temp]
# def fake_mess_rrho_partition_function(geo, freqs, hform0, temps):
#     """
#     """
#     def _pf_arrays(geo, freqs, temps):
#         temps_tuple = ()
#         lnq_tuple = ()
#         dlnqdt_tuple = ()
#         d2lnqdt2_tuple = ()
#         #####
#         q_tuple = ()
#         dqdt_tuple = ()
#         d2qdt2_tuple = ()
#         dlnqdlnt_tuple = ()
#         entropy = ()
#         enthalpy = ()
#         heat_cap = ()
#         for temp in temps:
#             temp_range = numpy.arange(temp-20, temp+50, .01)
#             lnq_temp_dct = rrho_partition_function(
#                 geo, freqs, temp_range, nlog=1)
#             lnq_func, dlnqdt_func, d2lnqdt2_func = pf_polys(
#                 lnq_temp_dct, order=3)
#             temps_tuple += (temp,)
#             lnq_tuple += (lnq_func(temp),)
#             dlnqdt_tuple += (dlnqdt_func(temp),)
#             d2lnqdt2_tuple += (d2lnqdt2_func(temp),)
#             #####
#             q_temp_dct = rrho_partition_function(
#                 geo, freqs, temp_range, nlog=0)
#             lnq_lnt_dct = rrho_partition_function(
#                 geo, freqs, temp_range, nlog=2)
#             q_func, dqdt_func, d2qdt2_func = pf_polys(
#                 q_temp_dct, order=3)
#             _, dlnqdlnt, _ = pf_polys(
#                 lnq_lnt_dct, order=2)
#             q_tuple += (q_func(temp),)
#             dqdt_tuple += (dqdt_func(temp),)
#             d2qdt2_tuple += (d2qdt2_func(temp),)
#             d2lnqdlnt2_tuple += (d2lnqdlnt2_func(temp),)
#
#             heat_cap_i = heat_capacity_from_pf(
#                 lnq_func, dlnqdlnt_func, d2lnqdt2_func, temp)
#             entropy_i = entropy_from_pf(q_func, dqdt_func, temp)
#             enthalpy_i = enthalpy_from_pf(q_func, dqdt_func, temp)
#             heat_cap += (heat_cap_i,)
#             entropy += (entropy_i,)
#             enthalpy += (enthalpy_i,)
#         print('Partition Function')
#         _inf1 = zip(temps_tuple, q_tuple, dqdt_tuple, d2qdt2_tuple)
#         _inf2 = zip(temps_tuple, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple)
#         _inf3 = zip(temps_tuple, entropy, enthalpy, heat_cap)
#         print('Partition Function')
#         for temp, y_val, dydx, d2ydx2 in _inf1:
#             print(temp, y_val, dydx, d2ydx2)
#         print('LOG Partition Function')
#         for temp, y_val, dydx, d2ydx2 in _inf2:
#             print(temp, y_val, dydx, d2ydx2)
#         print('Temp', 'Heat Cap', 'Entropy', 'Enthalpy')
#         for temp, heat_cap_i, entropy_i, enthalpy_i in _inf3:
#             print(temp, heat_cap_i, entropy_i, enthalpy_i)
#         return temps_tuple, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple
#
#     formula_str = automol.geom.formula_string(geo)
#     temps.append(298.15)
#     pf_arrays = _pf_arrays(geo, freqs, temps)
#     pf_str = mess_io.writer.pf_output(formula_str, *pf_arrays)
#     rundir = 'tmp'
#     print(pf_arrays)
#     # with tempfile.TemporaryDirectory() as rundir:
#     thermp_script_str = autorun.SCRIPT_DCT['thermp']
#     _, thermp_output_strs = thermp_direct(
#         thermp_script_str, rundir,
#         pf_str, formula_str, hform0, temps[:-1])
#     print(thermp_output_strs[0])
#     print(temps[:-1])
#     csh_t_dct = thermp_io.reader.properties_temp_dct(thermp_output_strs[0])
#     print(csh_t_dct)
#     # write_mess_output(formula_string, pf_arrays, rundir)


def from_ln_partition_function(lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple):
    """Translate partition functions as natural logs into
        nonlogged partition functions
    """
    pf_tuple = ()
    dqdt_tuple = ()
    d2qdt2_tuple = ()
    for lnq, dlnqdt, d2lnqdt2 in zip(lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple):
        pf_tuple += (numpy.exp(lnq),)
        dqdt_tuple += (numpy.exp(lnq) * dlnqdt,)
        d2qdt2_tuple += (numpy.exp(lnq) * (dlnqdt**2 + d2lnqdt2),)
    return (pf_tuple, dqdt_tuple, d2qdt2_tuple,)


def to_ln_partition_function(pf_tuple, dqdt_tuple, d2qdt2_tuple):
    """Translate partition functions to natural logs from
       nonlogged partition functions
    """
    lnq_tuple = ()
    dlnqdt_tuple = ()
    d2lnqdt2_tuple = ()
    for pf_val, dqdt, d2qdt2 in zip(pf_tuple, dqdt_tuple, d2qdt2_tuple):
        lnq_tuple += (numpy.log(pf_val),)
        dlnqdt_tuple += (dqdt / pf_val,)
        d2lnqdt2_tuple += (d2qdt2/pf_val - dqdt**2/pf_val**2,)
    return (lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple,)


def additive_pf_combination_at_temp(pf_arrays_lst, weight_lst, idx):
    """additively combine nonlog pfs, for one temperature idx, where
        each pf is given a weight
    """
    pf_tot_t = 0
    dqdt_tot_t = 0
    d2qdt2_tot_t = 0
    for weight_i, (pf_i, dqdt_i, d2qdt_i) in zip(weight_lst, pf_arrays_lst):
        pf_i_t, dqdt_i_t, d2qdt_i_t = pf_i[idx], dqdt_i[idx], d2qdt_i[idx]
        pf_tot_t += weight_i * pf_i_t
        dqdt_tot_t += weight_i * dqdt_i_t
        d2qdt2_tot_t += weight_i * d2qdt_i_t
    return (pf_tot_t, dqdt_tot_t, d2qdt2_tot_t,)


def weights_at_temp(pf_arrays_lst, hf_lst, temps, idx):
    """ Calculate the weight of each conformer given a list
        of partition functions and 0 K heats
        of formations for thtose conformers and a temperature
    """
    pfs_t_lst = []
    temp = temps[idx]
    hf_lst = [hf_i * phycon.EH2KJ * 1000 for hf_i in hf_lst]
    knt = phycon.KB * phycon.NAVO * temp
    for pf_array in pf_arrays_lst:
        pfs_t_lst.append(pf_array[0][idx])

    denominator = 0.
    h_zero = min(hf_lst)
    for q_val_i, h_val_i in zip(pfs_t_lst, hf_lst):
        exponent = numpy.exp((h_zero - h_val_i) / knt)
        denominator += q_val_i * exponent

    weight_lst = []
    for q_val_i, h_val_i in zip(pfs_t_lst, hf_lst):
        exponent = numpy.exp((h_zero - h_val_i) / knt)
        weight_lst.append(q_val_i * exponent / denominator)

    return weight_lst


def boltzmann_pf_combination(ln_pf_arrays_lst, hf_lst):
    """combine pfs
    """
    pf_arrays_lst = []
    for ln_pf_array in ln_pf_arrays_lst:
        temps, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple = ln_pf_array
        pf_arrays_lst.append(
            from_ln_partition_function(
                lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple))
    total_pf_arrays = ([], [], [])
    print('Weights:\n', 'Temperature (K)', 'Conformer Weight')
    for idx, temp in enumerate(temps):
        weight_lst = weights_at_temp(pf_arrays_lst, hf_lst, temps, idx)
        print(temp, '    ', '    '.join([f'{w:.3f}' for w in weight_lst]))
        pf_arrays_i = additive_pf_combination_at_temp(
            pf_arrays_lst, weight_lst, idx)
        total_pf_arrays[0].append(pf_arrays_i[0])
        total_pf_arrays[1].append(pf_arrays_i[1])
        total_pf_arrays[2].append(pf_arrays_i[2])
    final_ln_pf_arrays = to_ln_partition_function(*total_pf_arrays)
    final_ln_pf_arrays_with_temps = (temps, *final_ln_pf_arrays)

    return final_ln_pf_arrays_with_temps


def combine_pfs_additively(ln_pf_arrays_lst):
    """combine pfs additively
    """
    pf_arrays_lst = []
    for ln_pf_array in ln_pf_arrays_lst:
        temps, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple = ln_pf_array
        pf_arrays_lst.append(
            from_ln_partition_function(
                lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple))
    total_pf_arrays = ([], [], [])
    for idx, temp in enumerate(temps):
        weight_lst = [1]*len(pf_arrays_lst)
        pf_arrays_i = additive_pf_combination_at_temp(
            pf_arrays_lst, weight_lst, idx)
        total_pf_arrays[0].append(pf_arrays_i[0])
        total_pf_arrays[1].append(pf_arrays_i[1])
        total_pf_arrays[2].append(pf_arrays_i[2])
    final_ln_pf_arrays = to_ln_partition_function(*total_pf_arrays)
    final_ln_pf_arrays_with_temps = (temps, *final_ln_pf_arrays)

    return final_ln_pf_arrays_with_temps


def stereo_pf_combination(ln_pf_arrays_lst, hf_lst):
    """combine pfs
    """
    pf_arrays_lst = []
    for ln_pf_array in ln_pf_arrays_lst:
        temps, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple = ln_pf_array
        pf_arrays_lst.append(
            from_ln_partition_function(
                lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple))
    total_pf_arrays = ([], [], [])
    print('Stereo Weighting')
    print('Weights:\n', 'Temperature (K)', 'Conformer Weight')
    for idx, temp in enumerate(temps):
        weight_lst = stereo_weights(pf_arrays_lst, hf_lst, temps, idx)
        print(temp, '    ', '    '.join([f'{w:.3f}' for w in weight_lst]))
        pf_arrays_i = additive_pf_combination_at_temp(
            pf_arrays_lst, weight_lst, idx)
        total_pf_arrays[0].append(pf_arrays_i[0])
        total_pf_arrays[1].append(pf_arrays_i[1])
        total_pf_arrays[2].append(pf_arrays_i[2])
    final_ln_pf_arrays = to_ln_partition_function(*total_pf_arrays)
    final_ln_pf_arrays_with_temps = (temps, *final_ln_pf_arrays)

    return final_ln_pf_arrays_with_temps


def stereo_weights(pf_arrays_lst, hf_lst, temps, idx):
    """ Calculate the weight of each stereoisomer given a list
        of partition functions and 0 K heats
        of formations for thtose conformers and a temperature
    """
    pfs_t_lst = []
    temp = temps[idx]
    hf_lst = [hf_i * phycon.EH2KJ * 1000 for hf_i in hf_lst]
    knt = phycon.KB * phycon.NAVO * temp
    for pf_array in pf_arrays_lst:
        pfs_t_lst.append(pf_array[0][idx])

    weight_lst = []
    h_zero = min(hf_lst)
    for q_val_i, h_val_i in zip(pfs_t_lst, hf_lst):
        exponent = numpy.exp((h_zero - h_val_i) / knt)
        weight_lst.append(exponent)

    return weight_lst
