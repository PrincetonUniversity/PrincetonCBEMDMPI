#make_lj_test.py

import numpy as np

nbonds = 3

#make n1xn1xn1 grid of atoms
n1 = 5
natoms = n1**3
dim = 3
#set length of box
L = 10
#set diameter (not needed for LJ)
diameter = 1
#have map of types (0=A, 1=B)
type_lookup = ['A', 'B']
ntypes = len(type_lookup)

#create array of atom types
#make sure that A appears before B appears before C, etc,
#so that in indexing in LAMMPS and cbemd is the same
atom_type = np.random.randint(0, ntypes, size=natoms)
atom_type[:ntypes] = range(ntypes)

molecule_id = 1

#create grid of potisions
pos = np.zeros((natoms, dim))
ind = 0
for i in range(n1):
    for j in range(n1):
        for k in range(n1):
            pos[ind][0] = i / float(n1) * L
            pos[ind][1] = j / float(n1) * L
            pos[ind][2] = k / float(n1) * L
            ind += 1

#randomize velocities
vel = (np.random.rand(natoms, 3)-0.5)/2.0
#masses are all 1
mass = np.ones(ntypes)

#open and read to files
f1 = open('fene.xml','w')
f2 = open('fene.lammps','w')

f1.write('<?xml version="1.0" encoding="UTF-8"?> \n')
f1.write('<hoomd_xml version="1.4"> \n')
f1.write('<configuration time_step="0" dimensions="%d" natoms="%d" > \n' %(dim, natoms))
f1.write('<box lx="%d" ly="%d" lz="%d"/> \n' %(L, L, L))
f1.write('<position num="%d"> \n' %natoms)
for i in range(natoms):
    for j in range(dim):
        f1.write('%2.2f ' %pos[i][j])
    f1.write('\n')
f1.write('</position> \n')
f1.write('<velocity num="%d"> \n' %natoms)
for i in range(natoms):
    for j in range(dim):
        f1.write('%2.2f ' %vel[i][j])
    f1.write('\n')
f1.write('</velocity> \n')
f1.write('<mass num="%d"> \n' %natoms)
for i in range(natoms):
    f1.write('%2.2f ' %mass[atom_type[i]])
    f1.write('\n')
f1.write('</mass> \n')
f1.write('<diameter num="%d"> \n' %natoms)
for i in range(natoms):
    f1.write('%2.2f ' %diameter)
    f1.write('\n')
f1.write('</diameter> \n')
f1.write('<type num="%d"> \n' %natoms)
for i in range(natoms):
    f1.write(type_lookup[atom_type[i]])
    f1.write('\n')
f1.write('</type> \n')
f1.write('<bonds num="%d"> \n' %nbonds)
for i in range(nbonds):
    f1.write('bond1 %d %d \n' %(i, i+1))
f1.write('</bonds> \n')
f1.write('</configuration> \n')
f1.write('</hoomd_xml> \n')


f2.write('LAMMPS Description \n')
f2.write('\n')
f2.write('%d atoms \n' %natoms)
f2.write('%d bonds \n' %nbonds)
f2.write('0 angles \n')
f2.write('0 dihedrals \n')
f2.write('0 impropers \n')
f2.write('\n')

f2.write('%d atom types \n' %ntypes)
f2.write('1 bond types \n')
f2.write('\n')

f2.write('%2.2f %2.2f xlo xhi \n' %(0, L))
f2.write('%2.2f %2.2f ylo yhi \n' %(0, L))
f2.write('%2.2f %2.2f zlo zhi \n' %(0, L))
f2.write('\n')


f2.write('Atoms \n')
f2.write('\n')
for i in range(natoms):
    f2.write('%d ' %(i+1))
    f2.write('%d ' %molecule_id)
    f2.write('%d ' %(atom_type[i]+1))
    for j in range(dim):
        f2.write('%2.2f ' %pos[i][j])
    f2.write('\n')
f2.write('\n')

f2.write('Masses \n')
f2.write('\n')
for i in range(ntypes):
    f2.write('%d %2.2f' %(i+1, mass[i]))
    f2.write('\n')
f2.write('\n')

f2.write('Velocities \n')
f2.write('\n')
for i in range(natoms):
    f2.write('%d ' %(i+1))
    for j in range(dim):
        f2.write('%2.2f ' %vel[i][j])
    f2.write('\n')
f2.write('\n')

f2.write('Bonds \n')
f2.write('\n')
for i in range(nbonds):
    f2.write('%d ' %(i+1))
    f2.write('1 ')
    f2.write('%d %d' %(i+1, i+2))
    f2.write('\n')
f2.write('\n')

f1.close()
f2.close()


