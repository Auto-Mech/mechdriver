import autofile


f1 = 'Hfdb_0K.csv'
f2 = 'thermodb_0K.csv'

d1 = autofile.io_.read_file(f1).splitlines()[1:]
d2 = autofile.io_.read_file(f2).splitlines()[1:]

n1 = []

for line in d1:
    n1.append(line.split(',')[0])

for line in d2:
    name = line.split(',')[0]
    if name not in n1:
        if name + 'a' not in n1:
            print(line)
