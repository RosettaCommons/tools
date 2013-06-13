#!/usr/bin/python
#Author: David Masica
#January, 2004

"""
This script is designed to take two user input energy files, one from the ppk.pdb and one from the native.pdb and return the difference. It also gets rid of the annoying rossetta numbering system and replaces it with chain and pdb residue number
"""

import string
import sys

#Check that user input is correct

if len(sys.argv) != 3:
    print 'Usage: energy_diff.py <pdbfile1> <pdbfile2>'
    print 'Output: <energy_diff.pdbfile2.pdbfile1> = <pdbfile2> - <pdbfile1>'
    print 'Not enough arguments detected ....... exiting'
    sys.exit(1)

#This module strips off the \n characters and splits each line into list of strings

def fcn(mk_file):
    file = open(mk_file, 'r').readlines()

    stripfile = []
    for i in range(len(file)):
        stripfile.append(string.strip(file[i]))

    splitfile = []
    for i in range(len(stripfile)):
        splitfile.append(string.split(stripfile[i]))

    return splitfile

#This module removes all unescessary info from the input files

def get_table(energy_score):
    file = fcn(energy_score)
   
    new_file = []
    for i in range(len(file)):
        if len(file[i]) > 10:
            new_file.append(file[i])

    temp_file1 = []
    for i in range(len(new_file)):
        if new_file[i][0] != 'ATOM':
            temp_file1.append(new_file[i])

    temp_file2 = []
    for i in range(len(temp_file1)*0.5):
        temp_file2.append(temp_file1[i])

    energy_file = []
    for i in range(len(temp_file2)):
        if temp_file2[i][0] != 'res':
            energy_file.append(temp_file2[i])
    
    return energy_file

#This module takes input files and returns the difference and relables with pdb numbers and chains instead of rosseta numbers

def difference(pdbA, pdbB):
    pdb1 = get_table(pdbA)
    pdb2 = get_table(pdbB)
    file = fcn(pdbA)
        
    new_file = []
    for i in range(len(file)):
        if len(file[i]) > 10:
            new_file.append(file[i])

    temp_file1 = []
    for i in range(len(new_file)):
        if new_file[i][0] == 'ATOM' and new_file[i][2] == 'N':
            temp_file1.append(new_file[i])
    
    newfname = 'energy_diff' + '.' + sys.argv[2] + '.' + sys.argv[1]            
    file_out = open(newfname, 'w')
    file_out.write('ch res aa     Eatr   Erep   Esol   Eaa    Edun   Ehbnd  Epair  Eref   Eres' + '\n')
    for i in range(len(pdb1)):
        file_out.write(temp_file1[i][4] +'  '+ string.rjust(temp_file1[i][5], 2) +'  '+ pdb1[i][1] +'   '+ string.rjust(str(float(pdb2[i][2]) - float(pdb1[i][2])), 4) +'   '+ string.rjust(str(float(pdb2[i][3]) - float(pdb1[i][3])), 4) +'   '+ string.rjust(str(float(pdb2[i][4]) - float(pdb1[i][4])), 4) +'   '+ string.rjust(str(float(pdb2[i][5]) - float(pdb1[i][5])), 4) +'   '+ string.rjust(str(float(pdb2[i][6]) - float(pdb1[i][6])), 4) +'   '+ string.rjust(str(float(pdb2[i][7]) - float(pdb1[i][7])), 4) +'   '+ string.rjust(str(float(pdb2[i][8]) - float(pdb1[i][8])), 4) +'   '+ string.rjust(str(float(pdb2[i][9]) - float(pdb1[i][9])), 4) +'   '+ string.rjust(str(float(pdb2[i][10]) - float(pdb1[i][10])), 4) +'   '+ '\n' )
             
    
pdbA = sys.argv[1]
pdbB = sys.argv[2]
difference(pdbA, pdbB)


    
