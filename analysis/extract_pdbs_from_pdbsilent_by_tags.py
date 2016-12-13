#!/usr/bin/env python
#extract_pdbs_from_pdbsilent_by_tags.py

import sys

def read_tags( filename ) :
    file = open(filename, 'r')
    lines = []
    for line in file:
        lines.append(line[:-1])
    return lines
    
def print_a_pdb(silent_file, pdb_tag):
    ATOM = "ATOM"
    TER = "TER"
    print "making pdb for tag " + pdb_tag
    pdb_file = open(pdb_tag+".pdb", 'w')
    for line in silent_file:
        split_line = line.split()
        if( split_line[0] == ATOM):
            pdb_file.write(line)
        if( split_line[0] == TER):
            return
    

if __name__ == "__main__" :
    print '''usage: extract_pdbs_from_pdbsilent_by_tags.py silent_file tagsfile ; tagsfile is endline-delimited list of tags to extract; silent_file is a PDB-silent-file'''
    
    tagslist = read_tags(sys.argv[2])
    print "I am going to read these tags:"
    for each in tagslist:
        print each
        
    SEQUENCE = "SEQUENCE:"
    SCORE = "SCORE:"
    REMARK = "REMARK"
    ATOM = "ATOM"
    description = "description"
    TER = "TER"

    pdb_tag = ""
    
    silent_file = open(sys.argv[1])
    for line in silent_file:
        split_line = line.split()
        if (split_line[0] == SEQUENCE):
            #print line
            continue
        if (split_line[0] == SCORE):
            if (not (split_line[-1] == description)):
                pdb_tag = split_line[-1]
                #print line
                if pdb_tag in tagslist:
                    print_a_pdb(silent_file, pdb_tag)
                continue
            else:
                #print line
                continue
        if (split_line[0] == REMARK):
            #print line
            continue
        if (split_line[0] == ATOM):
            #print line
            continue
        if (split_line[0] == TER):
            #print line
            continue
        else:
            print "unhandled line"
            print line
            sys.exit()
            
