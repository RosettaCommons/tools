#!/usr/bin/python

import string
from os.path import exists,basename
from os import system, getcwd, chdir
from make_tag import make_tag, make_tag_with_dashes
from parse_options import get_resnum_chain
from get_sequence import get_sequences
from rna_server_conversions import get_all_stems, join_sequence
from sys import argv

user_input_runs = argv[1:]

info_file = 'favorites.txt'

names = []
sequence = {}
secstruct = {}
working_res = {}
native = {}
input_res = {}
extra_flags = {}
inpath = {}
lines = open( info_file ).readlines()

bps = ['au','ua','gc','cg','ug','gu']
fid_qsub = open( 'qsubMINI', 'w' )

for line in lines[ 1: ]:
    cols = string.split( line.replace( '\n','' ) )
    name = cols[0]
    if ( len( user_input_runs ) > 0 ) and (name not in user_input_runs): continue

    assert( not name in names ) # better be unique
    names.append( name )
    sequence [ name ]   = cols[1]
    secstruct[ name ]   = cols[2]
    working_res[ name ] = cols[3]
    native[ name ]      = cols[4]
    input_res[ name ]   = cols[5]
    extra_flags[ name ] = cols[6]
    inpath[ name ]      = cols[7]


    fasta_file = '%s/%s.fasta' % (inpath[name],name)
    fid = open( fasta_file, 'w' )
    sequences          = string.split( sequence[name], ',' )
    working_res_blocks = string.split( working_res[name], ',' )
    assert( len( sequences ) == len( working_res_blocks ) )
    for n in range( len( sequences ) ): fid.write( '>%s %s\n%s\n' % (name,working_res_blocks[n],sequences[n]) )
    fid.close()

    if native[ name ] != '-':
        starting_native = inpath[name]+'/'+native[name]
        assert( exists( starting_native ) )
        command = 'pdbslice.py %s -subset %s %s/%s ' % ( starting_native, string.join(working_res_blocks), inpath[name],name+'_' )
        system( command )
        working_native = inpath[name] + '/' + name + '_' + native[name]
        assert( exists( working_native ) )

        (sequences_native,all_chains,all_resnums) = get_sequences( working_native )
        assert( sequences_native ==  sequences )

    # create any helices.
    (sequence_joined, chainbreak_pos)           = join_sequence( sequence[name] )
    (secstruct_joined,chainbreak_pos_secstruct) = join_sequence( secstruct[name] )
    assert( chainbreak_pos == chainbreak_pos_secstruct )
    stems = get_all_stems( secstruct_joined, chainbreak_pos, sequence_joined  )

    resnums = []
    chains  = []
    for working_res_block in working_res_blocks: get_resnum_chain( working_res_block, resnums, chains )

    # helix creation block
    for i in range( len( stems ) ):
        stem = stems[i]
        helix_seq = ''; helix_resnum = []; helix_chains = []
        for bp in stem:
            helix_seq    += sequence_joined[ bp[0] - 1 ]
            helix_resnum.append( resnums[         bp[0] - 1 ] )
            helix_chains.append( chains[          bp[0] - 1 ]  )
        helix_seq += ' '
        for bp in stem[::-1]:
            helix_seq    += sequence_joined[ bp[1] - 1 ]
            helix_resnum.append( resnums[         bp[1] - 1 ] )
            helix_chains.append( chains[          bp[1] - 1 ] )
        helix_file =  '%s/%s_stem%d.pdb' % (inpath[name],name,(i+1))
        command = 'rna_helix.py -seq %s  -o %s -resnum %s' % ( helix_seq, helix_file, make_tag_with_dashes( helix_resnum, helix_chains) )
        print command
        system( command )
exit()

for name in names:
    dirname = name
    if not exists( dirname ): system( 'mkdir '+dirname )

    fid = open( '%s/README' % name, 'w' )
    fid.write( 'stepwise @flags -out:file:silent swm_rebuild.out\n' )
    fid.close()

    N = len( sequence )
    # tetraloop & dock site.
    input_res_start = [1,2,5,6,7,8,N-1,N]
    start_sequence = ''
    for m in input_res_start: start_sequence += sequence[ m-1 ]

    # figure out 'capping' helix.
    i = len( receptor_seq1[name] )
    j = 1
    input_res_helix = []
    helix_seq1 = ''
    helix_seq2 = ''
    print name
    while ( j <= len( start_helix_seq1[name] ) ):
        input_res_helix = [len(loop_seq[name]) + i] + \
                          input_res_helix + \
                          [len( loop_seq[name]) + len(receptor_seq1[name]) + j]
        assert( (receptor_seq1[name][i-1] + receptor_seq2[name][j-1]) in bps )
        helix_seq1 = receptor_seq1[name][i-1] + helix_seq1
        helix_seq2 = helix_seq2 + receptor_seq2[name][j-1]
        i -= 1
        j += 1
    assert( helix_seq1 == start_helix_seq1[name])
    assert( helix_seq2 == start_helix_seq2[name])
    input_res = input_res_start + input_res_helix
    helix_pdb = '%s_%s.pdb'  % (helix_seq1, helix_seq2)
    helix_pdb_with_path = '%s/%s' % (name, helix_pdb )
    if not exists( helix_pdb_with_path ):
        system( 'rna_helix.py -seq %s %s -o %s' % (helix_seq1,helix_seq2,helix_pdb_with_path ) )
    helix_pdb_sequence = get_sequence( helix_pdb_with_path )
    assert( helix_pdb_sequence == ( helix_seq1 + helix_seq2 ) )
    helix_sequence = ''
    for m in input_res_helix: helix_sequence += sequence[ m-1 ]
    assert( helix_sequence == helix_pdb_sequence )

    terminal_res = [1,len(loop_seq[name])] # from loop_res
    terminal_res.append( len(loop_seq[name]) + 1 )
    terminal_res.append( len(loop_seq[name]) + len(receptor_seq1[name] ) )
    terminal_res.append( len(loop_seq[name]) + len(receptor_seq1[name] ) + 1)
    terminal_res.append( N )

    fid = open( '%s/flags' % name, 'w' )
    system( 'cp %s %s/' % (start_pdb[name],name) )
    start_pdb_with_path = '%s/%s' % (name, basename(start_pdb[name]) )
    start_pdb_sequence = get_sequence( start_pdb_with_path )
    assert( start_pdb_sequence == start_sequence )

    fid.write( '-s %s %s\n' % (basename(start_pdb[name]), helix_pdb) )
    fid.write( '-fasta %s.fasta\n' % name )
    fid.write( '-terminal_res %s  \n' % make_tag( terminal_res ) )
    fid.write( '-extra_min_res %s \n' % make_tag( extra_min_res ) )
    fid.write( '-cycles 1000\n' )
    fid.write( '-nstruct 50\n' )
    # these should be read in from extra_benchmark_flags.txt
    fid.write( '-score:rna_torsion_potential RNA11_based_new\n' )
    fid.write( '-chemical::enlarge_H_lj\n' )

    if len( native_pdb[ name ] ) > 0:
        native_sequence = get_sequence( native_pdb[name] )
        assert( native_sequence == sequence )
        system( 'cp %s %s/' % (native_pdb[name],name) )
        fid.write( '-native %s\n' % basename( native_pdb[name] ) )

    if len( extra_flags[name] ) > 0 : fid.write( '%s\n' % extra_flags[name] )
    fid.close()

    CWD = getcwd()
    chdir( name )
    system( 'rosetta_submit.py README SWM 50' )
    chdir( CWD )

    fid_qsub.write( 'cd %s; source qsubMINI; cd %s\n' % ( name, CWD ) )

fid_qsub.close()

