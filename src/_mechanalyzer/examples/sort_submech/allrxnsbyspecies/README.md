Sort a mechanism according to the specified criteria running:
```
mechanalyzer sortmech
```
In this example, add input specifications:
```
mechanalyzer sortmech -m heptane_cut_mech.txt -s heptane_cut_species.csv
-i sort_singlespc.dat (or sort_singlespc_submech.dat) 
```

For additional info on input specifications,
```
mechanalyzer sortmech --help
```

Explanation of sorting:
sort\_singlespc.dat will extract all reactions involving C2H4
sort\_singlespc\_submech.dat will extract the submechanism involving C2H4,
i.e., also reactions with stoichiometry involving its decomposition and oxidation


