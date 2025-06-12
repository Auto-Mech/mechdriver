Umbrella
% if geo:
  Geometry[angstrom]   ${natom}
${geo}
% endif
  Group                ${group}
  Axis                 ${axis}
  ReferenceAtom        ${ref_atom} 
  Potential[kcal/mol]  ${npotential}
${potential} 
End  ! Umbrella
