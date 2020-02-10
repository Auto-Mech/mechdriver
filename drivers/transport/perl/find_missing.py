import sys

fname = sys.argv[1]

with open(fname,'r') as datfile:
  data = datfile.readlines()

for i in range(len(data)):
  if 'MISSING' in data[i]:
    print(data[i-6].strip()) 
