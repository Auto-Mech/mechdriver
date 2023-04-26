""" Install mechdriver
"""

from distutils.core import setup


setup(
    name='mechdriver',
    version='0.7.1',
    packages=[
        'drivers',
        'mechroutines',
        'mechroutines.es',
        'mechroutines.es._routines',
        'mechroutines.es.runner',
        'mechroutines.es.newts',
        'mechroutines.ktp',
        'mechroutines.models',
        'mechroutines.thermo',
        'mechroutines.proc',
        'mechroutines.trans',
        'mechroutines.trans._routines',
        'mechlib',
        'mechlib.amech_io',
        'mechlib.amech_io.parser',
        'mechlib.amech_io.printer',
        'mechlib.amech_io.reader',
        'mechlib.amech_io.writer',
        'mechlib.filesys',
        'mechlib.reaction'
    ],
    package_dir={
        'drivers': 'drivers',
        'mechroutines': 'mechroutines',
        'mechlib': 'mechlib'
    },
    scripts=['bin/automech.py'],
    }
)
