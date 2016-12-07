#!/usr/bin/env python

import string
from os.path import exists,basename


hetatm_map = { '5BU':'  U', ' MG':' MG', 'OMC':'  C', '5MC':'  C', 'CCC':'  C', ' DC':'  C', 'CBR':'  C', 'CBV':'  C', 'CB2':'  C', '2MG':'  G', 'H2U':'H2U', 'PSU':'PSU', '5MU':'  U', 'OMG':'  G', '7MG':'  G', '1MG':'  G', 'GTP':'  G', 'AMP':'  A', ' YG':'  G', '1MA':'  A', 'M2G':'  G', 'YYG':'  G', ' DG':'  G', 'G46':'  G', ' IC':' IC',' IG':' IG' }

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u',
              '  A': 'a', '  C': 'c', '  G': 'g', '  U': 'u',
              ' MG': 'Z[MG]',' IC':'c[ICY]',' IG':'g[IGU]',
              'ROS': 'Z[ROS]','HOH':'w[HOH]', 'H2U': 'X[H2U]',
              'PSU': 'X[PSU]', '5MU': 'X[5MU]', 'FME': 'X[FME]'
              }

from subprocess import Popen, PIPE
import os
grep = Popen( ["grep", "-r", "IO_STRING", "%s/main/database/chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_nonnatural/" % os.environ[ "ROSETTA" ] ], stdout=PIPE )
awk = Popen( ["awk", "{print $2}"], stdin=grep.stdout, stdout=PIPE )
grep.stdout.close()
tlcs, err = awk.communicate()
for tlc in tlcs.split('\n'):
    longer_names[tlc] = "X[%s]" % tlc

def get_sequences( pdbname, removechain = 0 ):

    netpdbname = pdbname
    assert( exists(netpdbname))

    lines = open(netpdbname,'r').readlines()

    oldresnum = '   '
    oldchain = ''
    chain = oldchain
    resnum = oldresnum
    count = 0;
    fasta_line = ''

    sequences = []
    all_chains = []
    all_resnums = []

    sequence = ''
    chains = []
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

        if ( line[0:3] == 'TER' or ( not chain == oldchain ) ) and len( sequence ) > 0:
            sequences.append( sequence )
            all_chains.append( chains )
            all_resnums.append( resnums )
            sequence = ''
            chains   = []
            resnums  = []
            old_resnum = ''

        if (not (resnum == oldresnum and chain == oldchain) ):
            count = count + 1
            longname = line_edit[17:20]
            if longer_names.has_key(longname):
                sequence +=  longer_names[longname]
            else:
                sequence +=  'X'
            resnums.append( int(resnum) )
            chains.append( chain )
            oldresnum = resnum
            oldchain = chain

    if len( sequence ) > 0:
        sequences.append( sequence )
        if ( chain == ' ' ): chain = ''
        all_chains.append( chains )
        all_resnums.append( resnums )

    return ( sequences, all_chains, all_resnums )

def get_sequence( pdbname, removechain = 0 ):
    ( sequences, chains, resnums ) = get_sequences( pdbname, removechain )
    return sequences[0]


if __name__=='__main__':

    import argparse

    parser = argparse.ArgumentParser(description='Get sequence from pdb.')
    parser.add_argument('pdbname', help='pdbfile to get sequence from')
    parser.add_argument('--removechain', action='store_true')
    args=parser.parse_args()

    ( sequences, all_chains, all_resnums ) = get_sequences( args.pdbname, removechain = args.removechain )
    print string.join(sequences, '')

