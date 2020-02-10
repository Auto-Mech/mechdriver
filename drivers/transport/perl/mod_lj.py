import sys

# Get name of file that contains LJ parameters
ljfilename = sys.argv[1]

# Set the N2-N2 LJ parameter values
n2_epsilon = 226.263856206897
n2_sigma = 3.40607149425287

# Set conversion factors
kB_cmK = 0.69503

# Grab all the lines from the file
with open(ljfilename,'r') as ljfile:
  data = ljfile.readlines()

# Find the line containing the list of all computed LJ parameters
endsecstart = 0
for i in range(len(data)):
  if 'End' in data[i] or 'MISSING' in data[i]:
    endsecstart = i + 2

print(endsecstart)
# Grab the info from the line in the LJ output, parse, and reorganize as necessary
newlines = []
for i in range(endsecstart,len(data)):
  tmp = data[i].strip().split()
  if len(tmp) > 4:
    name = tmp[0]
    epsilon = ( float(tmp[2])**2 / n2_epsilon ) / kB_cmK
    sigma = ( ( 2 * float(tmp[1]) ) - n2_sigma ) 
    newstr = '{0:15s}  2  {1:>4.3f}  {2:>4.3f}'.format(name,epsilon,sigma)
  else:
    name = tmp[0]
    newstr = '{0:15s}  2  x.xxx  x.xxx'.format(name)
  newlines.append(newstr)

# Set the text at the top of the file
#print('''! THEORETICAL TRANSPORT FOR SYNGAS MECHANISM
#!
#! Prepared by Moore, Jasper, and Klippenstein
#! 9/20/18
#! This file follows CHEMKIN's transport database format, with columns showing 
#! the species name followed by:
#! 1) A geometry parameter, indicating 
#!       0: an atom
#!       1: a linear molecule
#!       2: a nonlinear molecule
#! 2) epsilon (K), the Lennard-Jones well depth
#! 3) sigma (A), the Lennard-Jones collision diameter
#!
#! Parameters for A+N2 were calculated. The pure gas A+A parameters listed below were 
#! obtained by reversing the combining rules using the computed values for N2+N2.
#!
#! Best fit 12/6 LJ parameters to the "exact classical" diffusion coefficient in 
#! Jasper, Kamarchik, Miller, & Klippenstein, JCP, 141, 124313 (2014)
#H                0  541.672  1.530  0.000  0.666  0
#H2               1  304.690  2.190  0.000  0.775  280
#!
#! Calculated using the MP2/a'dz PES and the "1-D optimizations" method of 
#! Jasper & Miller, C&F 161, 101-110 (2014)''')


print('''! Species         Shape    LJ-depth  LJ-diam   DiplMom   Polzblty  RotRelaxNum Data     

! Name            Index    epsilon/k_B sigma     mu        alpha     Zrot      Source''')   

# Print the remaining lines containing organized information
for x in newlines:
  print(x)

