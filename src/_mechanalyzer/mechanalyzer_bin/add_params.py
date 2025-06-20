""" Add params to reaction lines for file to parse correctly
"""

import sys

# Set name of mechanism file
name = sys.argv[1]

# Read lines of file
with open(name, 'r') as f:
    flines = f.readlines()

# Generate new lines by altering original lines as needed
new_lines = []
for line in flines:
    sline = line.strip()
    if sline:
        if all(symb not in sline for symb in ('!', 'END')):
            new_lines.append('{0:<60s}{1:<40s}'.format(sline, '1.0  0.0  0.0'))
        else:
            new_lines.append(sline)

# Write out new lines
for line in new_lines:
    print(line)
