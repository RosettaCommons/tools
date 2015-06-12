#!/usr/bin/python

import string
from os.path import exists,basename
from parse_tag import parse_tag

longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u',
              '  A': 'a', '  C': 'c', '  G': 'g', '  U': 'u',
              ' MG': 'Z[MG]',' IC':'c[ICY]',' IG':'g[IGU]',
              }

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

        if line_edit[0:4] == 'ATOM':
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
    return sequences[0]

def get_sequences_for_res( pdbname, input_res, removechain = 0 ):
    sequences, chains, resnums = get_sequences( pdbname, removechain )
    input_resnums, input_chains = parse_tag(input_res)
    if all ( (c is None or not c.isalpha()) for c in input_chains ):
        input_chains = [c for chain in chains for c in chain]
    subsequences = []
    for sequence, subchains, subresnums in zip(sequences, chains, resnums):
        subsequence = ''
        for residue in zip(sequence, subchains, subresnums):
            if residue[2] not in input_resnums:
                continue
            if residue[1] not in input_chains:
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
    
    ( sequences, all_chains, all_resnums ) = get_sequences( args.pdbname, removechain = args.removechain )
    print string.join(sequences, '')
    
