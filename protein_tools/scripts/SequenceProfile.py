#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## Authors: Florian Richter
## Modified by Nick Randolph (Jan 2021)

################################################################################################
# Sequence profile for a list of structures vs a template
# To get the % of each residue at each site:
# SequenceProfile.py
# -l is the input list
# -t is the template to align to
# -c is the ID of chain of interest if multiple chains present
# NOTE: IF LOOKING FOR PROFILE OF MULTIPLE CHAINS, RESULTS MIGHT NOT BE ACCURATE! RUN SCRIPT 
#       MULTIPLE TIMES USING -c FOR EACH INDIVIDUAL CHAIN!
# python ~/install/perl/tools/SequenceProfile.py -l list.txt -t original_structure.pdb -c chainID
################################################################################################

import sys
import re


def AA3LetTo1Let(aa):

    if aa == 'ALA':
        return 'A'
    elif aa == 'CYS':
        return 'C'
    elif aa == 'ASP':
        return 'D'
    elif aa == 'GLU':
        return 'E'
    elif aa == 'PHE':
        return 'F'
    elif aa == 'GLY':
        return 'G'
    elif aa == 'HIS':
        return 'H'
    elif aa == 'ILE':
        return 'I'
    elif aa == 'LYS':
        return 'K'
    elif aa == 'LEU':
        return 'L'
    elif aa == 'MET':
        return 'M'
    elif aa == 'ASN':
        return 'N'
    elif aa == 'PRO':
        return 'P'
    elif aa == 'GLN':
        return 'Q'
    elif aa == 'ARG':
        return 'R'
    elif aa == 'SER':
        return 'S'
    elif aa == 'THR':
        return 'T'
    elif aa == 'VAL':
        return 'V'
    elif aa == 'TRP':
        return 'W'
    elif aa == 'TYR':
        return 'Y'
    else:
        return 'unknown'






def get_coordinates(struct,ChainID):

    struct = struct.replace("\n","")
    struct_file = open(struct,'r')
    struct_lines = struct_file.readlines()
    struct_file.close()
    atom_reading_flag = 0

    all_coordinates = {}
    cur_coordinates = {}
    cur_res = -50000
    for line in struct_lines:

        if (atom_reading_flag == 0) and line[0:4] == 'ATOM':
            atom_reading_flag = 1
        if (atom_reading_flag == 1 ) and line[0:4] != 'ATOM':
            atom_reading_flag = 0
            all_coordinates[cur_res] = cur_coordinates

        if atom_reading_flag:
            cols = line.split()
            if ( cols[4] == ChainID ) or ( ChainID == ''):
                if cur_res != int(cols[5]):
                    if cur_res != -50000:
                        all_coordinates[cur_res] = cur_coordinates
                    cur_res = int(cols[5])
                    cur_coordinates = {}
                    cur_coordinates['type'] = cols[3]
                    cur_coordinates['chain'] = cols[4]
                    cur_coordinates['active'] = 1  #keep track of whether we count this residue
                cur_coordinates[cols[2]] = (float(cols[6]), float(cols[7]), float(cols[8]))
               #     elif (not ca_only) and line[13:14] != 'H':
               #         cur_coordinates[cols[2]] = (float(line[31:38]), float(line[39:46]), float(line[47:54]))

    return all_coordinates


class SequenceProfile:

    def __init__(self, res_input):
        self.wt_pos = []
        self.res_id = []
	self.num_sequences = 0
        self.mutations = {}

        for res in res_input:
            self.wt_pos.append( AA3LetTo1Let(res_input[res]['type']) )
	    self.res_id.append( res )


    def add_struct( self, res_coords ):

        if len( res_coords ) != len( self.wt_pos ):
            print("Error: input structure and wt structure doesn't have the same size")
            print("res coords:", len(res_coords), "wt_pos", len(self.wt_pos))
            sys.exit()
        i = -1
        self.num_sequences = self.num_sequences + 1
        for res in res_coords:

            i = i + 1
            this_res = AA3LetTo1Let(res_coords[res]['type'])

            if this_res != self.wt_pos[ i ]:

                if not self.mutations.has_key( i ):
                    self.mutations[ i ] = {}
                if not self.mutations[ i ].has_key( this_res ):
                    self.mutations[i][this_res] = 0

                self.mutations[ i ][this_res] = self.mutations[i][this_res] + 1


    def get_outstring( self ):

        outstring = ""
        firstMut = 0
	for i in range( len( self.wt_pos ) ):
            if self.mutations.has_key( i ):
                if firstMut == 0:
		    outstring = "RosettaRes (pdbRes): % Mut\n"
		    firstMut = firstMut + 1
		
		outstring = outstring + self.wt_pos[i] + str(i+1) + " (" + self.wt_pos[i] + str(self.res_id[i]) + "): "
                for res in self.mutations[ i ]:

                    freq = float(self.mutations[ i ][res]) / float( self.num_sequences )
                    #outstring = outstring + str( freq ) + " " + res + ",  "
                    outstring = outstring + "%.2f " % freq + res + ",  "

                outstring = outstring + "\n"

        if outstring == "":
            outstring = "no mutations found"

        return outstring



FileList = []
Listfile = ''
template = ''
outfile = ""
Singlefile = ''
listmode = 0
ChainID = ''

CommandArgs = sys.argv[1:]

for arg in CommandArgs:
    if arg == '-t':
        template = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-l':
        Listfile = CommandArgs[CommandArgs.index(arg)+1]
        listmode = 1
    elif arg == '-out':
        outfile = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-s':
        Singlefile = CommandArgs[CommandArgs.index(arg)+1]
    elif arg == '-c':
	ChainID = CommandArgs[CommandArgs.index(arg)+1]



if ( (not Listfile and not Singlefile) or (not template) ):
     print("""
Usage:  python ~/install/perl/tools/SequenceProfile.py -l list.txt -t original_structure.pdb -c chainID

list.txt is a list of pdb files to be compared (e.g. ls *_DE_1.pdb > list.txt)
original_structure.pdb is the sequence against which all others are compared, often the original pdb, or the input to design
chainID is the chain of interest for the profile

Note that the script cannot profile  multiple chains. In those cases you will want to rerun the script multiple times, changing the chainID for each designed chain.
     """)
     sys.exit()
elif(listmode):
    inlist = open(Listfile,'r')
    FileList = inlist.readlines()
    inlist.close()
    if outfile == ' ':
        outfile = Listfile + '.ana'
    print("Checking structures for %s structures in %s to template %s" % (len(FileList), Listfile, template))

else:
    FileList.append(Singlefile)

template_coords = {}

if template != '':
    template_coords = get_coordinates(template,ChainID)

seq_prof = SequenceProfile( template_coords )

outstring = ""
pymutstring = "select mutations, resi "

for struct in FileList:

    if len(struct)>1:                             # make sure this is not an empty line
        struct_coords = get_coordinates(struct,ChainID)

        seq_prof.add_struct( struct_coords )

        #print mutstring

outstring = seq_prof.get_outstring()
    #print outstring


if outfile == "":
    print(outstring)
    #print pymutstring

else:
    outf = open(outfile,'w')
    outf.write(outstring)
    outf.close()
