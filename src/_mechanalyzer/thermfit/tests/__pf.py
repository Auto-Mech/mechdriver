""" test thermfit.pf
"""
import numpy
import thermfit.pf


# [temps, logq, dq_dt, d2q_dt2]
PFS = (
    # Partition Function Set 1
    ((500.0, 67.50302, 0.040505, 4.23039e-05),
     (1000.0, 75.20304, 0.0354286, 5.36839e-04),
     (1500.0, 84.2003, 0.0673282, 7.793224e-04),
     (2000.0, 103.40404, 0.0937334, 9.837278e-04)),
    # Partition Function Set 2
    ((500.0, 65.9811, 0.010705, 1.21039e-06),
     (1000.0, 72.3157, 0.0158886, 1.70879e-05),
     (1500.0, 82.8571, 0.0271782, 2.79324e-05),
     (2000.0, 100.375, 0.0437624, 3.83778e-05))
)
PFA = (
    (500.0, 67.50302, 0.040505, 4.23039e-05),
    (1000.0, 75.20304, 0.0354286, 5.36839e-04),
    (1500.0, 84.2003, 0.0673282, 7.793224e-04),
    (2000.0, 103.40404, 0.0937334, 9.837278e-04))
PFA = tuple(zip(*PFA))
PFB = (
    (500.0, 65.9811, 0.010705, 1.21039e-06),
    (1000.0, 72.3157, 0.0158886, 1.70879e-05),
    (1500.0, 82.8571, 0.0271782, 2.79324e-05),
    (2000.0, 100.375, 0.0437624, 3.83778e-05))
PFB = tuple(zip(*PFB))
PFS = (PFA, PFB)
COEFFS = (0.5, 1.2)
OPERATORS = ('multiply', 'divide')
COMBO = (74.625572, 0.0315206, 0.00043288878)
WEIGHTS = (0.8, 0.2)
HFS = (0., .2)
WEIGHTS_COMP = (0.9484, 0.0516)


def test__combine():
    """ test thermfit.pf.combine
    """

    combined_pfs = thermfit.pf.combine(PFS, COEFFS, OPERATORS)
    print(combined_pfs)


def test_to_and_from_ln_partition_function():
    """ test thermfit.pf.combine
    """

    _, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple = PFS[0]
    nonlog_pf_array = thermfit.pf.from_ln_partition_function(
        lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple)
    log_pf_array = thermfit.pf.to_ln_partition_function(*nonlog_pf_array)
    assert numpy.allclose(PFS[0][1], log_pf_array[0])
    assert numpy.allclose(PFS[0][2], log_pf_array[1])
    assert numpy.allclose(PFS[0][3], log_pf_array[2])


def test__additive_combo():
    pf_arrays_lst = (PFS[0][1:], PFS[1][1:])
    assert numpy.allclose(
        thermfit.pf.additive_pf_combination_at_temp(pf_arrays_lst, WEIGHTS, 1),
        COMBO)


def test__compute_weights():
    pf_array_lst = []
    for pf_i in PFS:
        temps, lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple = pf_i
        pf_array_lst.append(thermfit.pf.from_ln_partition_function(
            lnq_tuple, dlnqdt_tuple, d2lnqdt2_tuple))
    assert numpy.allclose(
        thermfit.pf.weights_at_temp(pf_array_lst, HFS, temps, 1),
        WEIGHTS_COMP)


def test__boltzman_partition_function():
    print(PFS)
    print(thermfit.pf.boltzmann_pf_combination(PFS, HFS))


if __name__ == '__main__':
    # test__combine()
    # test_to_and_from_ln_partition_function()
    # test__additive_combo()
    # test__compute_weights()
    test__boltzman_partition_function()
