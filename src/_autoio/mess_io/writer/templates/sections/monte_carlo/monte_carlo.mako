  MonteCarlo
    MoleculeSpecification   ${atom_list}
${flux_mode_str}\
    SymmetryFactor               ${sym_factor}
    DataFile                     ${data_file_name}
    ElectronicLevels[1/cm]       ${nlevels}
${levels}
% if ground_ene is not None:
    GroundEnergy[kcal/mol]       ${ground_ene}
% endif
% if reference_ene is not None:
    ReferenceEnergy[kcal/mol]    ${reference_ene}
% endif
## % if nfreqs > 0:
##   NonFluxionalFrequencies[1/cm]  ${nfreqs}
## % endif
% if ref_config_file_name:
  NoHessian
  ReferenceConfiguration           ${ref_config_file_name} 
% endif
% if excluded_volume_factor is not None:
  ExcludedVolumeFactor             ${excluded_volume_factor}
% endif
% if use_cm_shift:
  UseCMShift
% endif
