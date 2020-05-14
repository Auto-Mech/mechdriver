''' Install moldriver
'''
from distutils.core import setup


setup(name='moldriver',
      version='0.1.0',
      packages=['drivers',
                'routines',
                'routines.es',
                'routines.es._routines',
                'routines.es.runner',
                'routines.pf',
                'routines.pf.ktp',
                'routines.pf.thermo',
                'routines.pf.messf',
                'routines.pf.runner',
                'lib',
                'lib.amech_io',
                'lib.amech_io.parser',
                'lib.amech_io.writer',
                'lib.filesys',
                'lib.phydat',
                'lib.reaction'])
