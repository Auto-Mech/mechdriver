import os
#import numpy as np
import ioformat
from elstruct import reader


PATH = os.path.dirname(os.path.realpath(__file__))
DAT_PATH = os.path.join(PATH, 'data')

TEST_ZMA_G09 = (
    ('H', (None, None, None), (None, None, None), (None, None, None)), 
    ('X', (0, None, None), ('R1', None, None), 
     (1.8897261254578281, None, None)), 
    ('H', (0, 1, None), ('R2', 'A2', None), 
     (1.8285783674902638, 1.3532876358248318, None)), 
    ('N', (0, 1, 2), ('R3', 'A3', 'D3'), 
     (2.340558087208302, 1.5707963267948966, 2.836836401935811)), 
    ('H', (3, 0, 1), ('R4', 'A4', 'D4'), 
     (1.9390668745935322, 1.769221535363778, 0.0)), 
    ('H', (3, 0, 4), ('R5', 'A5', 'D5'), 
     (1.9390574259629048, 1.7695229362723048, 4.41228200660534))
)

TEST_ZMA_MOLPRO2015 = (
    ('H', (None, None, None), (None, None, None), (None, None, None)),
    ('X', (0, None, None), ('R1', None, None), 
     (1.8897261254578281, None, None)), 
    ('H', (0, 1, None), ('R2', 'A2', None), 
     (1.714926458852979, 1.3532881594236073, None)), 
    ('N', (0, 1, 2), ('R3', 'A3', 'D3'), 
     (2.4339672495896827, 1.6144347940825108, 2.8368372222405593)), 
    ('H', (3, 0, 1), ('R4', 'A4', 'D4'), 
     (1.9415046212953728, 1.7543438775808802, -0.10358703643511545)), 
    ('H', (3, 0, 4), ('R5', 'A5', 'D5'), 
     (1.9415046212953728, 1.7543630762026523, 4.426790282905346))
)


def test_g09():
    """ Tests the Gaussian09 reader (same as Gaussian16, I think)
    """
        
    fname = 'g09_zma.inp'
    inp_str = ioformat.pathtools.read_file(DAT_PATH, fname)
    prog = 'gaussian09'
    zma = reader._reader.inp_zmatrix(prog, inp_str)
    assert zma == TEST_ZMA_G09  # might fail due to numerical precision...


def test_molpro2015():
    """ Tests the Molpro2015 reader (same as other versions)
    """ 
    
    fname = 'molpro2015_zma.inp'
    inp_str = ioformat.pathtools.read_file(DAT_PATH, fname)
    prog = 'molpro2015'
    zma = reader._reader.inp_zmatrix(prog, inp_str)
    assert zma == TEST_ZMA_MOLPRO2015  # might fail due to numerical precision...


if __name__ == "__main__":
    test_g09()
    test_molpro2015()
    
