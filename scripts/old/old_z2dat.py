# import the path where the script is executed
import os    # to find pathway os
import re    # to use the re set (finds numbers in strings)
import array as arr # to use numeric arrays
cwd = os.getcwd()
#print(cwd)
# get the name of the file and open it
filename = "x2z_out.dat"
#print(filename)
path = '%s/%s'%(cwd,filename)
#path = '/home/chimica2/lpratali/Python/aim/to_convert.txt'
#print(path)

# find whether the molecule is linear
file = open(path,'r')
first_line = file.readline()
file.close()
if first_line.find("is nonlinear") != -1:
    ilin = 0
elif first_line.find("is plane") != -1:
    ilin = 0
elif first_line.find("is linear") != -1:
    ilin = 1
# step 1: find the symmetry number and save it
substring = "rotational symmetry number"
linecount = 0
line_symmetry = []
with open(path,'r') as myfile:
    for line in myfile:
        linecount += 1
        if line.find(substring) != -1:
            line_symmetry.append(line.rstrip('\n'))    # rstrip skippa di scrivere \n alla fine

string_elements = line_symmetry[0].split()    # split the string in different elements
symm_number = 'default'
for i in range(len(string_elements)):
#    print(string_elements[i].isdigit())
    if string_elements[i].isdigit() != False:
        symm_number = '%s'%(string_elements[i])

# find the line Z-Matrix and the line R1  = to set the boundaries where to read the zmatrix
# or just keep the if loop until you find a blank line
linecount = 0
dihedrals = []
with open(path,'r') as myfile:
    for line in myfile:
        linecount += 1
        # find lines delimiting z matrix
        if line.find('Z-Matrix:') != -1:
            line_in = linecount
#            print(line_in)
        if line.find('R1  =') != -1:
            line_fin = linecount-2
#            print(line_fin)
        # find rhe dihedrals
        if line.find('Rotational bond dihedral angles') != -1:
#            print(line)
            string_dih = re.split(' |,',line.rstrip('\n'))
#            print(string_dih)
            for x in string_dih:
                if x.find('D') != -1:
                    dihedrals.append(x)
#print(dihedrals)
str_dih = ''
for x in dihedrals:
    str_dih = '%s  %s'%(str_dih,x)
N_hind = len(dihedrals)
nosmp = (N_hind*N_hind)*5+1
N_atoms = line_fin-line_in
# select the N of sampling points 
#controllo che le linee contengano quello che devono
#file = open(path,'r')
#lines = file.readlines()
#print(lines[line_in],end='')
#print(lines[line_fin],end='')


elements = ['C','H','N','O','B','S','P','F','Cl','X']
indices = arr.array('i',[0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
atomic_N = arr.array('i',[6, 1, 7, 8, 5, 16, 15, 9, 17, 0])
linecount = -1
species = []
newline = []
final_intcoor = []
name_coor = []
R_coor = []
A_coor = []
D_coor = []
value_coor = []
k = 0
line_elements = []
dummy_coor = []
dummy_value = []
# add the indices to the species indices
with open(path,'r') as myfile:
    for line in myfile:
        linecount += 1
        if linecount >= line_in and linecount < line_fin:
            line_elements = [x.strip() for x in line.rstrip('\n').split(',')]
            element = '%s'%line_elements[0]
            dummy_flag = 0
            if element == 'X':
                # save the coordinates
                if len(line_elements) > 1:
                    dummy_coor.append(line_elements[2])
                if len(line_elements) > 3:
                    dummy_coor.append(line_elements[4])
                if len(line_elements) > 5:
                    dummy_coor.append(line_elements[6])
                
        if linecount > line_fin and linecount < line_fin+N_atoms:
            coordinates = [x.strip() for x in line.rstrip('\n').replace('=',' ').split()]
            kk = 0
            for x in coordinates:
                if ( kk % 2 ) == 0:
                    flag_dummy = 0
                    for z in dummy_coor:
                        if z == x :
                            flag_dummy = 1
                    if flag_dummy == 0:
                        if len(x) < 3:
                            name_coor.append('%s '%x)
                        else:
                            name_coor.append(x)
                else:
                    if flag_dummy == 0:
                        value_coor.append(x)
                    else:
                        dummy_value.append(x)
                kk += 1
# put the dot 
kk = 0
for x in value_coor:
    if value_coor[kk].find(".") == -1:
        value_coor[kk] = '%s.'%(value_coor[kk])
    kk += 1
kk = 0
for x in dummy_value:
    if dummy_value[kk].find(".") == -1:
        dummy_value[kk] = '%s.'%(dummy_value[kk])
    kk += 1
# reconstruct zmatrix
kk = 0                    
linecount = -1
bb = 0
with open(path,'r') as myfile:
    for line in myfile:
        linecount += 1
        if linecount >= line_in and linecount < line_fin:
            line_elements = [x.strip() for x in line.rstrip('\n').split(',')]
            # identify the species and overwrite the values
            element = '%s'%line_elements[0]
            for i in range(len(elements)):    
                if element == elements[i]:
                    indices[i] += 1
                    species.append('%s%i'%(element,indices[i]))
                    line_elements[0] = (species[len(species)-1])
            # if you have also the second element: check second and third element of the vector
            if len(line_elements) > 1:
                conn_1 = int(line_elements[1])
                line_elements[1] = species[conn_1-1]
            if len(line_elements) > 3:
                conn_2 = int(line_elements[3])
                line_elements[3] = species[conn_2-1]
                # add a space after the angle if there are fewer than 3 char.
            if len(line_elements) > 5:
                conn_3 = int(line_elements[5])
                line_elements[5] = species[conn_3-1]
            # save the line
            # put extra spaces 
            kk = 0
            for x in line_elements:
                if len(x) < 3:
                    line_elements[kk] = '%s '%line_elements[kk]
                kk += 1
            #correct dummy atoms
            if element == 'X':
                if len(line_elements) > 1:
                    line_elements[2] = dummy_value[bb]
                    while len(line_elements[2]) < 3 :
                        line_elements[2] = '%s '%(line_elements[2])
                    bb += 1
                if len(line_elements) > 3:
                    line_elements[4] = dummy_value[bb]
                    while len(line_elements[4]) < 3 :
                        line_elements[4] = '%s '%(line_elements[4])
                    bb += 1
                if len(line_elements) > 5:
                    line_elements[6] = dummy_value[bb]
                    bb += 1

            text_line = '  '.join(line_elements)
            full_line = '%s\n'%(text_line)
            newline.append(full_line)

kk = 0
for x in name_coor:
    if name_coor[kk].find("R") != -1:
        value_coor[kk] = str(float(value_coor[kk])*0.529177)
        R_coor.append('%s   %s\n'%(name_coor[kk],value_coor[kk]))
    elif name_coor[kk].find("A") != -1:
        A_coor.append('%s   %s\n'%(name_coor[kk],value_coor[kk]))
    elif name_coor[kk].find("D") != -1:
        check = 0
        for x in dihedrals:
            if name_coor[kk].find(x) != -1:
                check = 1
        if check != 1:
            D_coor.append('%s   %s\n'%(name_coor[kk],value_coor[kk]))
    kk +=1

kk = 0
for x in name_coor:
    check = 0
    for x in dihedrals:
        if name_coor[kk].find(x) != -1:
            check = 1
    if check == 1:
        D_coor.append('%s   %s\n'%(name_coor[kk],value_coor[kk]))
    kk += 1
final_intcoor = R_coor + A_coor + D_coor

# count the atoms
kk = 0
N_elec = 0
for x in indices:
    N_elec = N_elec + (x*atomic_N[kk])
    kk += 1
#print(N_elec)
if (N_elec % 2) == 0:
    spin = 1
else:
    spin = 2
# tot n of atoms: minus dummy
n_dummy = indices[-1]
N_atoms_net = N_atoms - n_dummy

# write file
path_name = ('%s/out_name.tmp'%cwd)
name_new = open(path_name,'r')
filename_new = name_new.read().rstrip('\n')
#print(filename_new)
name_new.close()
path_new = '%s/%s'%(cwd,filename_new)
file_new = open(path_new,'w')
file_new.write('nosmp dthresh ethresh\n')
file_new.write('%i  1.0  0.00001\n'%(nosmp))
file_new.write('\n')
file_new.write('ntau number of sampled coordinates\n')
file_new.write('%i\n'%(N_hind))
file_new.write('-->taumn,taumx sampling interval\n')
for x in dihedrals:
    if len(x) > 2:
        file_new.write('%s 0. 360.\n'%(x))
    else:
        file_new.write('%s  0. 360.\n'%(x))
file_new.write('\n')
file_new.write('nhind\n')
file_new.write('%i\n'%(N_hind))
file_new.write('-->namehind,hindmn,hindmx,nhindsteps\n')
for x in dihedrals:
    if len(x) > 2:
            file_new.write('%s 0. 360. 12 1\n'%(x))
    else:
        file_new.write('%s  0. 360. 12 1\n'%(x))
file_new.write('\n')
file_new.write('natom natomt ilin\n')
file_new.write('%i %i %i\n'%(N_atoms_net,N_atoms,ilin))
file_new.write('\n')
file_new.write('charge  spin  atomlabel\n')
file_new.write('0  %i \n'%(spin))
#file_new.write('0  %i\n'%(spin))
for x in newline:
    file_new.write(x)
file_new.write('\n')
file_new.write('intcoor\n')
for x in final_intcoor:
    file_new.write(x)
file_new.write('\n')
file_new.write('SymmetryFactor\n')
file_new.write('%s.\n'%symm_number)
file_new.write('\n')
file_new.write('nelec\n')
file_new.write('1\n')
file_new.write('0. %i.\n'%(spin))
file_new.write('\n')
file_new.write('end')
file_new.write('\n')

file_new.close()

