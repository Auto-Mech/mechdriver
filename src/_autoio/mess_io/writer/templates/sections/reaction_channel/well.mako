Well ${well_label} 
## Well Extension Cap Section
% if well_cap is not None:
  WellExtensionCap[kcal/mol]    ${well_cap}
% endif
  Species
${well_data}\
## Zero Energy Section
% if zero_ene is not None:
      ZeroEnergy[kcal/mol]      ${zero_ene}
% endif
  End  ! Species
## Energy Transfer Sections
% if edown_str is not None:
${edown_str}
% endif
% if collid_freq_str is not None:
${collid_freq_str}
% endif
End  ! Well\
