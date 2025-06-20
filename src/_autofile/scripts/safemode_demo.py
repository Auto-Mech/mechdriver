""" demonstrates use of the safemode setting
"""
import autofile
from autofile import fs

autofile.turn_off_safemode()

PFX = '/lcrc/project/PACC/AutoMech/data/save/'

print(autofile.safemode_is_on())

for zma_fs in fs.iterate_managers(
        PFX, ['REACTION', 'THEORY', 'TRANSITION STATE', 'CONFORMER'],
        'ZMATRIX'):
    print(zma_fs[0].path())
