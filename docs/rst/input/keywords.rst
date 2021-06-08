
Keywords
--------

Electronic Structure Driver Task Block
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Specifies electronic structure tasks

The order of the tasks in the input block corresponds to the order the tasks
are conducted in the order they are calculated.

Any number of tasks can be given in the file.


Explanation of keywords

init: geom
find: ts
conf: pucker, samp, energy, grad, hess, vpt2, prop, opt
hr: scan, grad, hess, energy, vpt2, reopt
tau: samp, energy, grad, hess
rpath: scan, energy, grad, hess:

| init: task where no info possibly exists using InChI or input geometry
| find: TS searching algorithm
| conf: conformer of species/TS defined torsional angles of internal rotors and ring torsions
| hr: all geometries along the motion of an internal rotor

| energy: calculates electronic energy
| grad: calculates molecular gradient in Cartesian coordinates
| hess: calculates molecular hessian in Cartesian coordinates
| vpt2: calculates VPT2 equations and saves several components
| prop: calculates molecular properties including dipole moment and polarizability
| opt: optimizes some geometry
| reopt: optimizes the minimum on the HR potential if a new min is found
| samp: optimizes conformers according Monte Carlo sampling routine
| scan: scans across some specified internal coordinate

| key descriptions:
| inplvl = electronic structure method where geometry is read
|          either for a guess (searching for scan) or for a calc on top (energy calc)
| runlvl = electronic structure method to calculate desired data for
| tors_model = model used to construct constraints for scan

for electronic structure method, these should correspond to those that defines blocks in 
theory.dat


kTPDriver and ThermDriver Task Blocks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The order given in the task block is not used anywhere. We use our own
internal order: write -> run -> fit

| write_mess: Reads and processes all information according to models specified
|             and then builds the required MESS file
|             give path to where the MESS file is saved in run-filesys
| run_mess: Run the MESS file existing in the filesystem
| run_fits: Fit the rate constants (thermochem) to functional forms and output
|           mechanism files with this data

kin_model: conditions model (kin mod) to use that is defined in models.dat
spc_model: conditions model (spc mod) to use that is defined in models.dat
overwrite: overwrite the MESS file currently in the file system
nprocs: number of processers to use during the MESS calculation
inpname:

Transport Driver Tasks
~~~~~~~~~~~~~~~~~~~~~~

| onedmin: Run OneDMin to calculate LJ epsilon and sigma
| write_transport: Build a ChemKin tranpsport file

Process Driver Tasks
~~~~~~~~~~~~~~~~~~~~

| freqs:
| energy: 
| geo:
| zmatrix:
| enthalpy:
| coeffs: 


species
~~~~~~~

inchi: InChI string (up to stereochemical layer [ref])
smiles: SMILES string (any form)

mult: spin multiplicity

charge: electric charge

inchikey: InChI key

sensitivity: sensitivity value for mechanism

tors_names: list of dihedral coords that define internal rotors

sym_factor: overall symmetry factor (external+internal)

hind_inc: increment of rotor potential (in degrees)

mc_nsamp: number of samples for Monte Carlo sampling of conformers [give formula]

tau_nsamp: number of samples for Monte Carlo sampling of conformers [give formula]

etrans_nsamp: samples for OneDMin

bath: bath gas molecule species interacts with OneDMin

smin: minimum intermolecular distance for OneDMin

smax: maximum intermolecular distance for OneDMin

lj: lennard-jones params [used for mess writing]

edown: alpha, n for energy down model (mess writing)

active: active space variables 

zma_idx: index for zma in save filesystem

rxn_dirn: direction to look for TS in 

kt_pst: k(T) for PST theory

temp_pst: temperature for k(T)

n_pst: n parameter for PST model potential

pst_params: [pre-exponential factor, pre-exponential power] Assuming a PST model potential

Maybe include  examples for input?

theory
~~~~~~

program: electronic structure program to run

method: electronic structure method

basis: basis set

orb_res: orb restiction label

memory: memory per core, in GB

nprocs: number of cores to use

econv: energy convergence threshhold

gconv: geometry convergence threshhold

optcoords: coordinates to use for geometry optimization

three combos of orb_res labels: 'RR', 'UU', 'RU' 
where the first (second) label in each pair corresponds closed-shell (open-shell) species, where R andU refers to use of an restricted or unrestricted reference, respectively


program (if version does not exist, just try other version, need to see if parsing makes filesys okay)

maybe provide a link to elstruct for program (method/basis) availability


models
~~~~~~

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

