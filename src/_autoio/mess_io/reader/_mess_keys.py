""" Defines keywords used to parse MESS input files
"""

OVERALL = [
'Model'
]

MODEL = [
'EnergyRelaxationFlags',
'CollisionFrequency',
'EnergyRelaxation',
'BufferFraction',
'Well',
'Barrier',
'Bimolecular',
'OutputTemperatureStep[K]',
'OutputTemperatureMin[K]',
'OutputTemperatureSize',
'OutputReferenceEnergy[kcal/mol]',
'ThermodynamicDataOutput',
'RelativeTemperatureIncrement',
'LumpingScheme',
'WellSeparator'
]

WELL = [
'EnergyRelaxation',
'Species',
'Escape',
'CollisionFrequency',
'WellExtensionCap[kcal/mol]',
]

SPECIES = [
'RRHO',
'Union',
'Variational',
'Read',
'Atom',
'MonteCarlo',
'MonteCarloWithDummyAtoms',
]

BIMOLEC = [
'GroundEnergy[1/cm]',
'GroundEnergy[kcal/mol]',
'GroundEnergy[au]',
'GroundEnergy[kJ/mol]',
'Fragment',
'Dummy',
]

ATOMIC = [
'Mass[amu]',
'Name',
'ElectronicLevels[1/cm]',
'ElectronicLevels[kcal/mol]',
'ElectronicLevels[kJ/mol]',
'ElectronicLevels[eV]',
'ElectronicLevels[au]',
]

RRHO = [
'ElectronicEnergy[1/cm]',
'ElectronicEnergy[kcal/mol]',
'ElectronicEnergy[kJ/mol]',
'ElectronicEnergy[eV]',
'ElectronicEnergy[au]',
'ZeroEnergy[1/cm]',
'ZeroEnergy[kcal/mol]',
'ZeroEnergy[kJ/mol]',
'ZeroEnergy[eV]',
'ZeroEnergy[au]',
'ElectronicLevels[1/cm]',
'ElectronicLevels[kcal/mol]',
'ElectronicLevels[kJ/mol]',
'ElectronicLevels[eV]',
'ElectronicLevels[au]',
'InterpolationEnergyMax[kcal/mol]',
'InterpolationEnergyStep[1/cm]',
'ExtrapolationStep',
'Frequencies[1/cm]',
'FrequencyScalingFactor',
'FrequencyDegeneracies',
'AreFrequenciesHarmonic',
'Anharmonicities[1/cm]',
'Rotor',
'HinderedRotorBundle',
'Core',
'Tunneling',
'InfraredIntensities[km/mol]',
'GraphPerturbationTheory',
#'SymmetryFactor',  # this conflicts with Core and is rarely used
]

