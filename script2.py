""" temporary script
"""
# input arguments
SCRIPT_STR = "#!/usr/bin/env bash\npsi4 >> stdout.log &> stderr.log"
PROG = 'psi4'
METHOD = 'rhf-mp2'
BASIS = 'sto-3g'
ZMA = (('C', (None, None, None), (None, None, None)),
       ('O', (0, None, None), (5.0556787897, None, None)),
       ('H', (0, 1, None), (3.902303346, 1.91162422312, None)),
       ('H', (0, 1, 2), (3.902303346, 1.91162422312, 2.10849736274)),
       ('H', (0, 1, 2), (3.90149076409, 1.90209472540, 4.1958413349)),
       ('H', (1, 0, 2), (3.472333961, 1.86909054925, 5.22893662580)))
VAR_DCT = {
    (1, 0): 'R1', (2, 0): 'R2', (3, 0): 'R3', (4, 0): 'R4',
    (5, 1): 'R5', (2, 0, 1): 'A2', (3, 0, 1): 'A3', (4, 0, 1): 'A4',
    (5, 1, 0): 'A5', (3, 0, 1, 2): 'D3', (4, 0, 1, 2): 'D4',
    (5, 1, 0, 2): 'D5'
}
MULT = 1
CHARGE = 0

# run some things
# 0. D5 is our active coordinate (could have multiple)
# 1. determine ranges for each coordinate based on periodicity (determine in
#    automol); allow manual specification
# 2. [input] the number of monte carlo samples
# 3. pick random numbers in the range for each active coordinate
# 4. for each combination set the values and get the (optimized,
#    constrained-optimized, or frozen) energy
# 5. return the unique geometries with their energies
#    a. find a good way to check uniqueness of geometries
# NOTE:
