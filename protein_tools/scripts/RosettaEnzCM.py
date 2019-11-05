import os
import sys
import Bio
import argparse
from Bio import AlignIO
from Bio import SeqIO
from Bio.SeqUtils import seq3
from Bio.SeqUtils import seq1
import pandas as pd
import pyrosetta
pyrosetta.init()

# Authors: Dr. Jason S. Fell, Dr. Timothy Coulther, Dr. Stephanie C. Contreras,
# Dr. Steve J. Bertolani and Dr. Justin B. Siegel 
#
# See: Bertolani SJ, Siegel JB (2019) A new benchmark illustrates that integration of geometric constraints 
# inferred from enzyme reaction chemistry can increase enzyme active site modeling accuracy. 
# PLoS ONE 14(4): e0214126. https://doi.org/10.1371/journal.pone.0214126
#
# This script identifies the catalytic residues in target sequences and generate distance constraints for 
# CA-CA, CB-CB, CA-CB, and CB-CA atoms between the catalytic residues. 
#
# This script calls Biopython and PyRosetta.
#
# To run this script: python RosettaEnzCM.py -c CATALYTIC-RESIDUE-INFORMATION -a ALIGNED-FASTA.fasta -n NAME 
#
# Within your working directory  place these files:
#    ** aligned fasta file with all template and target
#    ** all template pdbs
#    ** a file with template names (no .pdb) and catalytic residues in three-letter-code+position comma-separated:
#                IE:
#                >>  TEMPLATE1 RES1,RES2,RES3
#                >>  TEMPLATE2 RES1,RES2,RES3
# 
# Ensure all names are consistent in each file and respective filenames.
#
# The first part of this script will determine the catalytic residues for each target sequence and will save 
# the catalytic residue information for each template in a .data file. If your template sequences are mis-aligned
# the script will stop working.
#
# The second part of this script will then calculate the atomic distance constraints for each template 
# and save distance constraints for each target in TARGET.dist_csts file.
#
# NOTE: THIS SCRIPT WILL FAIL IF A TARGET RESIDUE IS DETERMINED TO BE ALANINE OR GLYCINE!
#
# This script saves many .txt, .csv and .data files for the user to use/call later.

##############PART i: Catalytic Residue Identifier#############################################

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--cat', dest='catres', type=str, help="File holding catalytic residue information.", required = True)
parser.add_argument('-a', '--aligned', dest='aligned', type=str, help='Aligned fasta file.', required = True)
parser.add_argument('-n', '--name', dest='name', type=str, default = 'RosettaEnzCM', help='Name of output files (optiponal)', required = False)
args = parser.parse_args()

print('Catalytic Residue File = %s' % (args.catres))

data = {}
for fasta in SeqIO.parse("%s" % (args.aligned), "fasta"): # Counts each residue position
    seq_list = []
    counter = 0
    for char in fasta.seq:
        if char == '-':
            seq_list.append(char)
        else:
            counter += 1
            seq_list.append('%s%s' % (char, counter))
    data[fasta.id] = seq_list

df1 = pd.DataFrame(data) #A newdataframe with numbered residues with residue one-letter-codes
temp_res = {}

with open('%s' % (args.catres), 'r') as f: # Read in catalytic information about templates. #
    for line in f:
        list2 = []
        a = line.strip()
        z = a.split(" ")
        with open("%s.data" %(z[0]), 'w') as q: # Writes individual catalytic data to file for each template. #
            q.write(z[1])
        q.close()
        b = a.replace(",", " ")
        c = b.split(" ")
        list = c # Strips spaces and commas from input file. #
        for i in range(1,len(list)):
            residue = str(list[i])
            one = seq1(residue[0:3])
            position = str(residue[3:])
            one_res = str("%s%s" % (one, position))
            list2.append(one_res)
        temp_res[list[0]] = list2 # Generates dictionary with template and one-letter-code residue info. #

template = []
positions = []
temp_pos = {}

#Confirmation that catalytic residues within the templates are aligned within the dataframe#

df1.to_csv("%s_info.csv" % (args.name))
for i in temp_res:
    print("Template = " + i)
    template.append(i)
    posit = [] # Temporary Positions List
    for j in temp_res[i]:
        print("Residue = " + j)
        posit.append(df1[df1[i] == j].iloc[0].name)
        if df1[df1[i] == j].iloc[0].name in positions:
            pass
        else:
            positions.append(df1[df1[i] == j].iloc[0].name)
    temp_pos[i] = posit

nomatchs = 0 # No Matches Counter
for i in temp_pos: # Prints message to confirm that templates match each other! #
    if temp_pos[i] == positions:
        print(i + " Matches")
    else:
        print(i + " Does Not Match!!!!!!")
        nomatchs += 1 # Adds to No matches counter

if nomatchs != 0:
    print('There was a mismatch! Recheck your alignment or your residues.')
    sys.exit()
        
print(positions)
res_count = len(positions)

family = {}

for i in df1: # Parsing family dataframe for each catalytic residue for each target. #
    name = i
    locations = []
    for j in positions:
        locations.append(df1[name].loc[j])
    family[name] = locations # Generates a new family dictionary of catalytic residues. #

A_G_count = 0 # Alanine and Glycine Counter

with open('%s_residues.txt' % (args.name), 'w') as g: # Writing target catalytic residue information in three-letter-code. #
    g.write('ID ')
    for n in range(1,res_count):
        g.write('Res%s Pos%s ' % (n, n))
    g.write('Res%s Pos%s\n' % (res_count, res_count))
    for i in family:
        if i in template:
            pass
        else:
            g.write(i)
            for j in family[i]:
                residue = j
                three = seq3(residue[0:1]) # Converts one- to three-letter-code. #
                position = str(residue[1:])
                three_res = str("%s %s" % (three, position))
                g.write(" %s" % (three_res.upper()))
                if residue[0:1] == 'A':
                    A_G_count += 1
                    print('Alanine Detected!')
                if residue[0:1] == 'G':
                    A_G_count += 1
                    print('Glycine Detected!')
            g.write('\n')
g.close()

if A_G_count != 0:
    print('Alanine/Glycine was detected! Recheck your alignment and templates.')
    sys.exit()
else:
    print('No Alanine or Glycine determined as catalytic residue!')

#########################PART ii: Catalytic constraints generator#############################

residue1 = []
residue2 = []
pairing = []
CA1CA2 = []
CB1CB2 = []
CA1CB2 = []
CB1CA2 = []
FullTemp = []

for i in os.listdir('.'): # Parsing thru data files. #
    if i.endswith('.data'):
        template_n = i[:-5]
        res_list = []
        with open(i, 'r') as f:
            for line in f:
                a = line.strip()
                b = a.replace(",", " ")
                c = b.split(" ")
                list = c
                for j in range(len(list)):
                    res = str(list[j])
                    position = int(res[3:])
                    res_list.append(position)
        pose = pyrosetta.rosetta.core.import_pose.pose_from_file('%s.pdb' % (template_n))    
        for k in range(0,len(res_list)): # Now we will begin calculating atomic distances. #
            pairing.append(k+1)
            if k <= len(res_list)-2:
                x = res_list[k]
                y = res_list[k+1]
            else:
                x = res_list[0]
                y = res_list[len(res_list)-1]    
            CA1 = pose.residue(x).xyz('CA')
            CB1 = pose.residue(x).xyz('CB')
            CA2 = pose.residue(y).xyz('CA')
            CB2 = pose.residue(y).xyz('CB')
            CA1CA2_xyz = CA1 - CA2
            CB1CB2_xyz = CB1 - CB2
            CA1CB2_xyz = CA1 - CB2
            CB1CA2_xyz = CB1 - CA2
            residue1.append(x)
            residue2.append(y)
            CA1CA2.append(CA1CA2_xyz.norm())
            CB1CB2.append(CB1CB2_xyz.norm())
            CA1CB2.append(CA1CB2_xyz.norm())
            CB1CA2.append(CB1CA2_xyz.norm())
            FullTemp.append(template_n)    #adding templates up
            
df = pd.DataFrame(zip(FullTemp,pairing,residue1,residue2,CA1CA2,CB1CB2,CA1CB2,CB1CA2), 
                  columns=['template','pairing','residue1','residue2','CA1CA2','CB1CB2','CA1CB2','CB1CA2'])

pairs = []
for i in set(pairing):
    pairs.append(i)
  
CACA_average = []
CBCB_average = []
CACB_average = []
CBCA_average = []

with open('%s_distances.txt' % (args.name),'w') as distance: # Calculates distances and saves to external file
    distance.write('Pair CACA CBCB CACB CBCA\n')
    for i in pairs:
        distance.write('%s ' % (i))
        df_copy = df[df['pairing'] == i]
        CACA_average.append(df_copy['CA1CA2'].mean())
        distance.write('%s ' % (df_copy['CA1CA2'].mean()))
        CBCB_average.append(df_copy['CB1CB2'].mean())
        distance.write('%s ' % (df_copy['CB1CB2'].mean()))
        CACB_average.append(df_copy['CA1CB2'].mean())
        distance.write('%s ' % (df_copy['CA1CB2'].mean()))
        CBCA_average.append(df_copy['CB1CA2'].mean())
        distance.write('%s\n' % (df_copy['CB1CA2'].mean()))
distance.close()       
listdata = []

d = pd.read_csv('%s_residues.txt' % (args.name), sep='\s+') # Opens last file from first step
listdata.append(d) # Creates new dataframe with catalytic identifiers
df_fam = pd.concat(listdata)

for i in df_fam['ID']: # Slices new DF for each target. #
    if i in template:
        pass
    else:
        print('Catalytic constraints for target: ' + i)
        with open('%s.dist_csts' % (i), 'w') as h: # Writes catalytic constraints to file. #
            df_famtest = df_fam[df_fam['ID']==i]
            for j in pairs:
                if j <= pairs[-2]: # Generates residue pairs based upon number of catalytic residues. #
                    res1 = df_famtest['Pos%s' % (j)].iloc[0]
                    res2 = df_famtest['Pos%s' % (j+1)].iloc[0]
                else: 
                    res1 = df_famtest['Pos%s' % (pairs[0])].iloc[0]
                    res2 = df_famtest['Pos%s' % (pairs[-1])].iloc[0]
                ca_ca_d = CACA_average[pairs.index(j)]
                cb_cb_d = CBCB_average[pairs.index(j)]
                ca_cb_d = CACB_average[pairs.index(j)]
                cb_ca_d = CBCA_average[pairs.index(j)]
                h.write('AtomPair CA %s CA %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(ca_ca_d, 2)))
                h.write('AtomPair CB %s CB %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(cb_cb_d, 2)))
                h.write('AtomPair CA %s CB %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(ca_cb_d, 2)))
                h.write('AtomPair CB %s CA %s SCALARWEIGHTEDFUNC 1000 HARMONIC %s 1.0\n' % (res1, res2, round(cb_ca_d, 2)))
        h.close


