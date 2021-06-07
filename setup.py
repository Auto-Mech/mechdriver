''' Install moldriver
'''
from distutils.core import setup


setup(name='moldriver',
      version='0.1.0',
      packages=['drivers',
                'mechroutines',
                'mechroutines.es',
                'mechroutines.es._routines',
                'mechroutines.es.runner',
                'mechroutines.pf',
                'mechroutines.pf.ktp',
                'mechroutines.pf.thermo',
                'mechroutines.pf.messf',
                'mechroutines.pf.runner',
                'mechlib',
                'mechlib.amech_io',
                'mechlib.amech_io.parser',
                'mechlib.amech_io.printer',
                'mechlib.amech_io.runner',
                'mechlib.amech_io.writer',
                'mechlib.filesys',
                'mechlib.reaction'])
