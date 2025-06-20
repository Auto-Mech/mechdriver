!===================================================
!  GLOBAL KEYWORDS
!===================================================
% if temperatures != '':
TemperatureList[K]                     ${temperatures}
% else:
Temperature(step[K],size)              ${temp_step}    ${ntemps}
% endif:
RelativeTemperatureIncrement           ${rel_temp_inc}
% if float_type == 'quadruple':
AtomDistanceMin[angstrom]              ${atom_dist_min}
!
FloatType                              dd\
% else:
AtomDistanceMin[angstrom]              ${atom_dist_min}\
% endif
