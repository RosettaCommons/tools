#!/usr/bin/env python

from __future__ import print_function
import string
from os.path import exists,basename
from parse_tag import parse_tag


hetatm_map = { '5BU':'  U', ' MG':' MG', 'OMC':'  C', '5MC':'  C', 'CCC':'  C', ' DC':'  C', 'CBR':'  C', 'CBV':'  C', 'CB2':'  C', '2MG':'  G', 'H2U':'H2U', 'PSU':'PSU', '5MU':'  U', 'OMG':'  G', '7MG':'  G', '1MG':'  G', 'GTP':'  G', 'AMP':'  A', ' YG':'  G', '1MA':'  A', 'M2G':'  G', 'YYG':'  G', ' DG':'  G', 'G46':'  G', ' IC':' IC',' IG':' IG', 'ZMP':'ZMP', 'YYG':'  G', '2MG':'  G','H2U':'  U' }

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u',
              '  A': 'a', '  C': 'c', '  G': 'g', '  U': 'u',
              ' MG': 'Z[MG]',' IC':'c[ICY]',' IG':'g[IGU]',
              'ROS': 'Z[ROS]','HOH':'w[HOH]', 'H2U': 'X[H2U]',
              'PSU': 'X[PSU]', '5MU': 'X[5MU]', 'FME': 'X[FME]',
	      'U33': 'X[U33]', '  I': 'X[INO]', 'BRU': 'X[5BU]' }

from subprocess import Popen, PIPE
import os
grep = Popen( ["grep", "-r", "IO_STRING", "%s/main/database/chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_nonnatural/" % os.environ[ "ROSETTA" ] ], stdout=PIPE )
awk = Popen( ["awk", "{print $2}"], stdin=grep.stdout, stdout=PIPE )
grep.stdout.close()
tlcs, err = awk.communicate()
for tlc in tlcs.decode('iso-8859-15').split('\n'):
    longer_names[tlc] = "X[%s]" % tlc

def get_sequences( pdbname, removechain = 0 ):

    netpdbname = pdbname
    assert( exists(netpdbname))

    lines = open(netpdbname,'r').readlines()

    oldchain = ''
    oldsegid = "    "
    oldresnum = '   '
    chain = oldchain
    segid = oldsegid
    resnum = oldresnum
    count = 0;
    fasta_line = ''

    sequences = []
    all_chains = []
    all_segids = []
    all_resnums = []

    sequence = ''
    chains = []
    segids = []
    resnums = []

    for line in lines:
        if (len(line)<20): continue
        line_edit = line
        if (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
            line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
            if (line_edit[12:14] == 'SE'):
                line_edit = line_edit[0:12]+' S'+line_edit[14:]
            if len(line_edit)>75:
                if (line_edit[76:78] == 'SE'):
                    line_edit = line_edit[0:76]+' S'+line_edit[78:]
        elif (line[0:6] == 'HETATM') & ( line[17:20] in hetatm_map.keys()):
            line_edit = 'ATOM  '+line[6:17]+ hetatm_map[line[17:20]] + line[20:]

        if (line[0:6] == 'HETATM') & (line[17:20] in longer_names.keys() ):
            line_edit = 'ATOM  '+line[6:]

        if line_edit[0:4] == 'ATOM' or line_edit[0:6]=='HETATM':
            resnum = line_edit[22:26].replace( ' ', '' )
            chain = line_edit[21]
            segid = line_edit[72:76]
        if ( line[0:3] == 'TER' or ( not chain == oldchain ) or ( not segid == oldsegid ) ) and len( sequence ) > 0:
            sequences.append( sequence )
            all_chains.append( chains )
            all_resnums.append( resnums )
            all_segids.append( segids )
            sequence = ''
            chains   = []
            resnums  = []
            segids = []
            old_resnum = ''

        if (not (resnum == oldresnum and chain == oldchain )):#and segid == oldsegid ) ):
            count = count + 1
            longname = line_edit[17:20]
            if longname in longer_names.keys():
                sequence +=  longer_names[longname]
            else:
                sequence +=  'X'
            resnums.append( int(resnum) )
            chains.append( chain )
            segids.append( segid )
            oldresnum = resnum
            oldchain = chain
            oldsegid = segid

    if len( sequence ) > 0:
        sequences.append( sequence )
        if ( chain == ' ' ): chain = ''
        all_chains.append( chains )
        all_resnums.append( resnums )
        all_segids.append( segids )

    return ( sequences, all_chains, all_resnums, all_segids )

def get_sequence( pdbname, removechain = 0, join = False ):
    sequences, chains, resnums = get_sequences( pdbname, removechain )
    if join:
        return ','.join(sequences)
    return sequences[0]

def get_sequences_for_res( pdbname, input_res, removechain = 0 ):
    sequences, chains, resnums, segids = get_sequences( pdbname, removechain )
    input_resnums, input_chains, input_segids = parse_tag(input_res)
    if all ( (c is None or not c.isalpha()) for c in input_chains ):
        input_chains = [c for chain in chains for c in chain]
    subsequences = []
    for sequence, subchains, subresnums, subsegids in zip(sequences, chains, resnums, segids ):
        subsequence = ''
        for residue in zip(sequence, subchains, subresnums, subsegids):
            if residue[2] not in input_resnums:
                continue
            if residue[1] not in input_chains:
                continue
            if residue[3] not in input_segids:
                continue
            subsequence += residue[0]
        if len(subsequence):
            subsequences.append(subsequence)
    return subsequences

if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Get sequence from pdb.')
    parser.add_argument('pdbname', help='pdbfile to get sequence from')
    parser.add_argument('--removechain', action='store_true')
    args=parser.parse_args()

    ( sequences, all_chains, all_resnums, all_segids ) = get_sequences( args.pdbname, removechain = args.removechain )
    print(''.join(sequences))
