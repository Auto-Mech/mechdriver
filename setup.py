''' Install moldriver
'''
from distutils.core import setup


setup(name='mechdriver',
      version='0.1.1',
      packages=['drivers',
                'mechroutines',
                'mechroutines.es',
                'mechroutines.es._routines',
                'mechroutines.es.runner',
                'mechroutines.pf',
                'mechroutines.pf.ktp',
                'mechroutines.pf.models',
                'mechroutines.pf.thermo',
                'mechroutines.pf.runner',
                'mechroutines.output',
                'mechroutines.trans',
                'mechroutines.trans._routines',
                'mechlib',
                'mechlib.amech_io',
                'mechlib.amech_io.parser',
                'mechlib.amech_io.printer',
                'mechlib.amech_io.runner',
                'mechlib.amech_io.writer',
                'mechlib.filesys',
                'mechlib.reaction'])
