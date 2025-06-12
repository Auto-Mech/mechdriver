Bimolecular ${bimolec_label} 
!---------------------------------------------------
  Fragment ${spc1_label}
${spc1_data}\
% if not isatom1:
      ZeroEnergy[kcal/mol]    0.0
% endif
  End  ! Frag1
!---------------------------------------------------
  Fragment ${spc2_label}
${spc2_data}\
% if not isatom2:
      ZeroEnergy[kcal/mol]    0.0
% endif
  End  ! Frag2
!---------------------------------------------------
  GroundEnergy[kcal/mol]    ${ground_ene}
End  ! Bimol\
