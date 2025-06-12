Barrier ${rxn_label}
${ts_data}\
## Zero Energy Section
% if zero_ene is not None:
    ZeroEnergy[kcal/mol]      ${zero_ene}
% endif
## Tunnel Section
% if tunnel != '':
${tunnel}
  End
% endif
End  ! Barrier\
