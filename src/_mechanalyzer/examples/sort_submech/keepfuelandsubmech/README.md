Sort a mechanism according to the specified criteria running:
```
mechanalyzer sortmech
```
In this example, add input specifications:
```
mechanalyzer sortmech -m heptane_cut_mech.txt -s heptane_cut_species.csv
```

For additional info on input specifications,
```
mechanalyzer sortmech --help
```

Explanation of sorting:
submechanism of C2H4 extracted (i.e., oxidation and decomposition)
then, all reactions below a given stoichiometry are kept;
the rest are deleted.
for additional info on available options: check sort.dat

