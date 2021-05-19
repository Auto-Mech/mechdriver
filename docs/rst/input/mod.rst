
models.dat
==========

File used to describe how the partition functions will be constructed

| kin model: specify the conditions by which rates./thermo is calculated and fit
|
| pressures: pressures used to calculate k(T,P) values with the master equation
| rate_temps: temperatures used to calculate k(T,P) values with the master equation
| therm_temps: temperatures used to calculate partition functions used for
|     thermochemical parameter determination
| rate_fit: methodology used to fit the rate constants to functional forms
|     fit_method: functional form to fit rate constants to
|     pdep_temps: list of temperatures at which to assess pressure dependence
|     pdep_tol: % error tolerance for determining if reaction is pressure dep
|     pdep_pval: pressure value to get rate constants for if no pressure dependence
|     pdep_plow: low pressure to assess pressure difference to
|     pdep_phigh: high pressure to assess pressure difference
|     arrfit_dbltol: % error tolerance for
|     troefit_params: list of parameters to fit in Troe expression
| therm_fit: methodology used to fit thermochemical parameters to NASA polynomials
|     ref_scheme: CBH scheme used to calculate energies
|     ref_enes: reference energies used to calcualte delta H

| Spc model: specify the means for calculating the constituent components of the
|            partition functions used in rate and thermochemistry calcualtion.
| 
|         ene: relative energies
|         rot: rotational pf, geometry and
|         vib: vibrational pf, frequencies
|         tors: vibrational pf, internal rotor degrees of freedom
|         symm: symmetry number
|         ts: handle transition states
| 
|     mod: model used to calculate the
|     geolvl: electronic structure level to read geometry
|     enelvl: electronic structure level to read energy
|     vpt2lvl: electronic structure level to get from VPT2 calc
|     tunnel: tunneling model
|     sadpt: Transition state theory for saddle points
|     nobar: Transition state theory method used for barrierless reactions
|     wells: manner by which build entrance-/exit-channel wells for a TS

theory levels correspond to ones provided in the theory.dat file

for enes

Example::

    ene = [0.2, lvl_cc_df]

