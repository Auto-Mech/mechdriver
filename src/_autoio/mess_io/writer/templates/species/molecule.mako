RRHO
## Core Section
${core}
## Frequencies Section
% if anharm == '':
    % if nfreqs > 0:
  Frequencies[1/cm]         ${nfreqs}
${freqs}
    % endif
% endif
## Electronic Levels Section
  ElectronicLevels[1/cm]    ${nlevels}
${levels}
## Hindered Rotor Section
% if hind_rot != '':
${hind_rot}\
% endif 
## Infrared Intensities
% if nintens > 0:
  InfraredIntensities[km/mol]  ${nintens}
${intens}
% endif
## Various Keywords
% if freq_scale_factor is not None:
  FrequencyScalingFactor    ${freq_scale_factor} 
% endif
% if use_harmfreqs_key: 
  AreFrequenciesHarmonic
% endif
