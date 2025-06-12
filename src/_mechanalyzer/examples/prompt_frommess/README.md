Compute branching fractions of prompt products 
from radical-generation and radical decomposition PESs
previously assembled and run with MESS
and obtain final rate constants in CHEMKIN format

run using the command
```
mechanalyzer promptcalc
```
in this example, add input specifications as:
```
mechanalyzer promptcalc -o rate.out
```
For additional info on input specifications,
```
mechanalyzer promptcalc --help
```

output:
bf\_model\_reaction\_promptproduct : branching fraction of the prompt product of interest (T,P) from given reaction
rates\_prompt.txt : rate constants fitted in PLOG format


