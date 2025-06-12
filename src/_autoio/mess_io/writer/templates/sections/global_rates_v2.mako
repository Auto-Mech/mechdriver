!===================================================
!  GLOBAL KEYWORDS
!===================================================
TemperatureList[K]                     ${temperatures}
PressureList[atm]                      ${pressures}
! 
ReferenceTemperature[K]                ${ref_temperature}
ReferencePressure[atm]                 ${ref_pressure}
!
ModelEnergyLimit[kcal/mol]             ${model_ene_limit}
EnergyStepOverTemperature              ${ene_stepover_temp}
EnergyCutoffOverTemperature            ${ene_cutoff_temp}
ExcessEnergyOverTemperature            ${excess_ene_temp}
!
ChemicalTolerance                      ${chem_tol}
ChemicalThreshold                      ${chem_thresh}
WellProjectionThreshold                ${well_projection_thresh}
WellReductionThreshold                 ${well_reduction_thresh}
!
TimePropagationLimit                   ${time_propagation_limit}
TimePropagationStep                    ${time_propagation_step}
!
WellExtension                          ${well_extension}
% if ground_ene_shift_max is not None:
GroundEnergyShiftMax[kcal/mol]         ${ground_ene_shift_max}
% endif
% if ped_spc_str is not None:
!
PEDOutput                              ${ped_outname}
PEDSpecies                             ${ped_spc_str}
% endif
% if hot_ene_str is not None:
!
HotEnergies[kcal/mol]                  ${nhot}
${hot_ene_str}
% endif
% if micro_out_params is not None:
!
MicroRateOutput                        ${ke_outname}
MicroEnerMin[kcal/mol]                 ${micro_out_params[0]}
MicroEnerMax[kcal/mol]                 ${micro_out_params[1]}
MicroEnerStep[kcal/mol]                ${micro_out_params[2]}
% endif
!
AtomDistanceMin[angstrom]              0.68793
!
RateOutput                             ${ktp_outname}
% if float_type == 'quadruple':
!
FloatType                              dd
% endif
