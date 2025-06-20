memory,${memory},m

${wfn_guess}

GEOMETRY_HERE

if (iterations.ge.0) then
   ${method}
   molpro_energy = energy + ${inf_sep_energy}
else
   molpro_energy = 10.0
endif
---
