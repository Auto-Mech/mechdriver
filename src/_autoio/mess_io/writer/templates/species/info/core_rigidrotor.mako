  Geometry[angstrom]        ${natom}
${geo}
  Core RigidRotor
    SymmetryFactor          ${sym_factor}
% if interp_emax is not None:     
    ZeroPointEnergy[1/cm]             0.0
    InterpolationEnergyMax[kcal/mol]  ${interp_emax}
% endif
% if anharm != '':
    % if nfreqs > 0:
    Frequencies[1/cm]         ${nfreqs}
${freqs}
    % endif
    Anharmonicities[1/cm]
${anharm}
    ## Rovibrational Coupling Section
    % if rovib_coups != '':
    RovibrationalCouplings[1/cm]
${rovib_coups}
    % endif
    ## Rotational Distortion Section
    % if rot_dists != '':
    RotationalDistortion[1/cm]
${rot_dists}
    End
    % endif
% endif 
  End  ! Core
