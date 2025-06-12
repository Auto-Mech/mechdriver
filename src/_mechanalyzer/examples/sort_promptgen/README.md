Extract prompt reactions from a mechanism according to the specified criteria running:
```
mechanalyzer sortmech
```
in this example, add input specifications as:
```
mechanalyzer sortmech -m DMM.CKI -t DMM.THERM
```

sort.dat contains:
filtering criteria to apply to select prompt reactions to analyze (see file for additional comments)

outmech.dat contains:
all reactions involving the selected prompt species
Comments are added if the reaction is analyzed, i.e., if it is
a radical decomposition or radical generation reaction

pes\_groups.dat should be copied in the inp directory to run calculations
it lists PES numbers to be considered according to the filtering criteria imposed in sort.dat (if any)

rxn\_prompt\_dh.out:
list of analyzed radical-generation reactions with the respective DH, the estimated "hot" DH,
and the total DH (of radical generation + radical decomposition reaction) computed from thermo for reaction filtering.
the "YES" and "NO" indicate if the reaction should be kept for analysis.
The selected reactions are the same as those specified in pes\_groups.dat

For additional info on input specifications,
```
mechanalyzer sortmech --help
```
