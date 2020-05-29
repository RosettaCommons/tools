#!/usr/bin/env python
from os import system,popen
import string
import argparse

#takes a pdb file, and the remodel range, and fills nterm and cterm with the residue numbers and amino acid.
#nterm will be from start of pdb to remodel range i.e (1,S)..(range[0],N)
#cterm will be from end of remodel range to end of pdb i.e (range[1],F)..(120,G) if structure is 120 aa long.
def get_aa_sequence_from_pdb(pdb_file,remodel_range,nterm,cterm):
    pdb = map(string.split,open(pdb_file,'r').readlines())
    amino_acids = {"ALA":"A","GLY":"G","ILE":"I","LEU":"L","PRO":"P","VAL":"V","PHE":"F","TRP":"W","TYR":"Y","ASP":"D","GLU":"E","ARG":"R","HIS":"H","LYS":"K","SER":"S","THR":"T","CYS":"C","MET":"M","ASN":"N","GLN":"Q"}
    resnum = 0
    for i in range(0, len(pdb)):
        if (pdb[i][0] == "ATOM"):
            if (pdb[i][5] != resnum):
                if (int(pdb[i][5]) <= remodel_range[0]):
                    resnum = pdb[i][5]
                    aa_ss = (resnum,amino_acids[pdb[i][3]])
                    nterm.append(aa_ss)
                if (int(pdb[i][5]) >= remodel_range[1]):
                    resnum = pdb[i][5]
                    aa_ss = (resnum,amino_acids[pdb[i][3]])
                    cterm.append(aa_ss)

def write_Nterm(nterm,ss,filename):
    out = open(filename, 'a')
    for i in range(0,len(nterm)-1):
        out.write(str(nterm[i][0])+" "+str(nterm[i][1])+" .\n")
    lastline = len(nterm)-1
    out.write(str(nterm[lastline][0])+" "+str(nterm[lastline][1])+" "+ss+"\n")

def write_Cterm(cterm,ss,filename):
    out = open(filename, 'a')
    out.write(str(cterm[0][0])+" "+str(cterm[0][1])+" "+ss+"\n")
    for i in range(1,len(cterm)-1):
        out.write(str(cterm[i][0])+" "+str(cterm[i][1])+" .\n")
    lastline = len(cterm)-1
    out.write(str(cterm[lastline][0])+" "+str(cterm[lastline][1])+" .")

#takes the filename to output, along with segment_info in form tuple(length_of_segment,secondary_structure)
#writes extention of length length_of_segment to filename 
def write_extension(filename, segment_info):
    out = open(filename, 'a')
    stub="0 x "+segment_info[1] +"\n"
    count = 0
    while count < segment_info[0]:
        out.write(stub)
        count += 1


#takes nterm and cterm data, to print out base, unchanging information, as well as data from three segments i,j,and k that are of tuple(length_of_segment,secondary_structure)
def write_blueprint(nterm,cterm,seg1,seg2,seg3,filename):
    write_Nterm(nterm,seg1[1],filename)
    write_extension(filename,seg1)
    write_extension(filename,seg2)
    write_extension(filename,seg3)
    write_Cterm(cterm,seg3[1],filename)

#takes in the imput PDB file, to extract sequence from, the remodel range, the length of three remodel segments that can have different secondary structures, and the base name for the output file.
#this function will iterate through all the lengths of the remodel segments, and call write_blueprint() to write each iteration
def write_all_blueprint(input_file,remodel_range,first_seg,second_seg,third_seg,tag):
    nterm = []
    cterm = []
    get_aa_sequence_from_pdb(input_file,remodel_range,nterm,cterm)
    count1 = first_seg[0]
    while count1 <= first_seg[1]:
        count2  = second_seg[0]
        while count2 <= second_seg[1]:
            count3 = third_seg[0]
            while count3 <= third_seg[1]:
                filename = tag+'_'+str(count1)+'_'+str(count2)+'_'+str(count3)+'.bp'
                print (filename)
                seg1 = (count1,first_seg[2])
                seg2 = (count2,second_seg[2]) 
                seg3 = (count2,third_seg[2])
                write_blueprint(nterm,cterm,seg1,seg2,seg3,filename)
                count3 +=1
            count2 +=1
        count1 += 1



parser = argparse.ArgumentParser(description="Takes in a pdb and makes a remodel blueprint file. Will generate multiple bluprints if --multi is true(default). WARNING: multi functionality is hard coded in the script. I wouldn't advise using it to make a basic backbone file, since it will spit out your filename with '{$tag}_0_0_0.bp'") 
parser.add_argument('PDB', type=str, nargs=1, help="PDB file to be made into blueprint")
parser.add_argument('-m', '--multi', type=bool, default = True, nargs=1, help="Print out multiple blueprint files with varied remodel lengths. WARNING: these numbers are hard coded, copy this script and change for own use")
parser.add_argument('-o','--out', type=str, nargs=1, default=["remodel"], help="Name of output file")

args = parser.parse_args()
input_file=args.PDB[0] #name of PDB to extract the sequence from
tag = args.out[0] #name of the output pdb

########################
#     FILE OPTIONS     #
########################

#blueprint files will be the same until hits the start_remodel residue, then it will output l*m*n blueprint files, based on ranges for remodeling residues. i.e first_seg_range = (0,3) second_seg_range = (2,3) third_seg_range = (0,3) will spit out a total of 32 files.
remodel_range = (55,63) #choose the residues from input PDB to be remodeled. everything within this range will be replaces with the 3 segments below. residues at listed positions will not change residue, but will me remodeled.
first_seg_range = (0,3,"H") #least residues added, most residues added, secondary structure, choose "H" for helix "L" for loop or "E" for sheet
second_seg_range = (2,3,"L") #least residues added, most residues added, secondary structure, choose "H" for helix "L" for loop or "E" for sheet
third_seg_range = (0,3,"H") #least residues added, most residues added, secondary structure, choose "H" for helix "L" for loop or "E" for sheet

if (args.multi == False):
    first_seg_range = second_seg_range = third_seg_range = (0,0,".")

write_all_blueprint(input_file,remodel_range,first_seg_range,second_seg_range,third_seg_range,tag)
