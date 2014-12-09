#!/usr/bin/python

import string
import argparse
from os.path import exists,basename,dirname
from os import system, getcwd, chdir
from make_tag import make_tag_with_conventional_numbering, make_tag_from_list_of_int_ranges
from parse_options import get_resnum_chain
from get_sequence import get_sequences
from rna_server_conversions import get_all_stems, join_sequence
from sys import argv, exit

parser = argparse.ArgumentParser(description='Setup benchmark for stepwise monte carlo')
parser.add_argument("info_file",       help='text file with information, in same directory as input_files/ (e.g., "../favorites.txt")')
parser.add_argument("user_input_runs", nargs='*',help='specify particular cases to run (default: run all in info_file)' )
default_extra_flags_benchmark = 'extra_flags_benchmark.txt'
parser.add_argument('-extra_flags', default=default_extra_flags_benchmark, help='Filename of text file with extra_flags for all cases.')
parser.add_argument('-nhours', default='16', type=int, help='Number of hours to queue each job.')
parser.add_argument('-legacy', default=False, help='Additional flag for setting up SWA runs.')
args = parser.parse_args()

# read in benchmark information & create any missing input files.
lines = open( args.info_file ).readlines()
relpath = dirname( args.info_file )
if len( relpath) > 0: relpath += '/'

names = []
sequence = {}
secstruct = {}
working_res = {}
native = {}
input_res = {}
extra_flags = {}
inpath = {}
fasta = {}
resnums = {}
chains  = {}
helix_files = {}
working_native = {}
input_pdbs = {}
terminal_res = {}
extra_min_res = {}

bps = ['au','ua','gc','cg','ug','gu']
fid_qsub = open( 'qsubMINI', 'w' )

for line in lines[ 1: ]:
    if line[0] == '#': continue
    if len( line.replace(' ','') ) < 5: continue
    cols = string.split( line.replace( '\n','' ) )
    name = cols[0]
    print 'Reading in info for: ', name
    if ( len( args.user_input_runs ) > 0 ) and (name not in args.user_input_runs): continue

    assert( not name in names ) # better be unique
    names.append( name )
    sequence [ name ]   = cols[1]
    secstruct[ name ]   = cols[2]
    working_res[ name ] = cols[3]
    native[ name ]      = cols[4]
    input_res[ name ]   = cols[5]
    inpath[ name ]      = relpath + '/input_files/' + basename( args.info_file.replace('.txt','' ) )
    extra_flags[ name ] = string.join( cols[6:] )
        
    
    if extra_flags[ name ] == '-': extra_flags[ name ] = ''
    sequences          = string.split( sequence[name], ',' )
    working_res_blocks = string.split( working_res[name], ',' )

    # create fasta
    fasta[ name ] = '%s/%s.fasta' % (inpath[name],name)
    if not exists( fasta[ name ] ):
        fid = open( fasta[ name ], 'w' )
        assert( len( sequences ) == len( working_res_blocks ) )
        for n in range( len( sequences ) ): fid.write( '>%s %s\n%s\n' % (name,working_res_blocks[n],sequences[n]) )
        fid.close()

    # store information on 'conventional' residue numbers and chains.
    resnums[ name ] = []
    chains[ name ] = []
    for working_res_block in working_res_blocks: get_resnum_chain( working_res_block, resnums[ name ], chains[ name ] )

    # helper function for PDB processing
    def slice_out( inpath_dir, prefix, pdb, res_string ):
        starting_native = inpath_dir+'/'+pdb
        assert( exists( starting_native ) )
        slice_pdb = prefix + pdb
        if not exists( slice_pdb ):
            command = 'pdbslice.py %s -subset %s %s ' % ( starting_native, res_string, prefix )
            system( command )
        assert( exists( slice_pdb ) )
        return slice_pdb

    # working_native
    assert( native[ name ] != '-' ) # for now, require a native, since this is a benchmark.
    prefix = '%s/%s_' % ( inpath[name],name)
    working_native[ name ] = slice_out( inpath[ name ], prefix, native[ name ], string.join( working_res_blocks ) )
    assert( string.join(sequences,'') == string.join(get_sequences( working_native[name] )[0],'') )

    # create starting PDBs
    input_pdbs[ name ] = []
    input_resnums = []
    input_chains  = []
    if input_res[ name ] != '-':
        input_res_blocks = string.split( input_res[ name ], ';' )
        for m in range( len ( input_res_blocks ) ):
            prefix = '%s/%s_START%d_' % ( inpath[name],name,m+1)
            input_pdb = slice_out( inpath[ name ], prefix, native[ name ],input_res_blocks[m] )
            input_pdbs[ name ].append( input_pdb )
            get_resnum_chain( input_res_blocks[m], input_resnums, input_chains )
    
    def get_fullmodel_number( reschain, resnums, chains):
        for m in range( len( resnums ) ):
            if ( resnums[m] == reschain[0] ) and ( reschain[1] == '' or chains[m] == reschain[1] ): return m+1
        return 0
    
    input_resnum_fullmodel = map( lambda x: get_fullmodel_number(x,resnums[name],chains[name]), zip( input_resnums, input_chains ) )


    # create secstruct if not defined
    if secstruct[ name ] == '-': 
        secstruct[ name ] = ''
        res_ranges = working_res[ name ].split(',')
        for res_range in res_ranges:
            res_str = make_tag_from_list_of_int_ranges( [ res_range ] )
            secstruct[ name ] += string.join(['.' for x in res_str.split()],'')
            secstruct[ name ] += ','
        secstruct[ name ] = secstruct[ name ][:-1]
        print 'Secstruct for '+name+': '+secstruct[ name ]


    # create any helices.
    helix_files[ name ] = []
    (sequence_joined, chainbreak_pos)           = join_sequence( sequence[name] )
    (secstruct_joined,chainbreak_pos_secstruct) = join_sequence( secstruct[name] )
    assert( chainbreak_pos == chainbreak_pos_secstruct )
    stems = get_all_stems( secstruct_joined, chainbreak_pos, sequence_joined  )

    for i in range( len( stems ) ):
        helix_file =  '%s/%s_HELIX%d.pdb' % (inpath[name],name,(i+1))

        stem = stems[i]
        helix_seq = ''; helix_resnum = [];
        for bp in stem:
            helix_seq    += sequence_joined[ bp[0] - 1 ]
            helix_resnum.append( bp[0] )
        helix_seq += ' '
        for bp in stem[::-1]:
            helix_seq    += sequence_joined[ bp[1] - 1 ]
            helix_resnum.append( bp[1] )

        already_in_input_res = False
        for m in helix_resnum:
            if m in input_resnum_fullmodel: already_in_input_res = True
        if already_in_input_res: continue

        helix_files[ name ].append( helix_file )
        input_resnum_fullmodel += helix_resnum

        if exists( helix_file ): continue
        command = 'rna_helix.py -seq %s  -o %s -resnum %s' % ( helix_seq, helix_file, \
            make_tag_with_conventional_numbering( helix_resnum, resnums[ name ], chains[ name ] ) )
        print command
        system( command )

    L = len( sequence_joined )
    terminal_res[ name ] = []
    extra_min_res[ name ] = []
    for m in range( 1, L+1 ):
        if ( m not in input_resnum_fullmodel ): continue
        prev_moving = ( m - 1 not in input_resnum_fullmodel ) and ( m != 1 )
        next_moving = ( m + 1 not in input_resnum_fullmodel ) and ( m != L )
        right_before_chainbreak = ( m == L or m in chainbreak_pos )
        right_after_chainbreak  = ( m == 1 or m - 1 in chainbreak_pos )
        if ( ( right_after_chainbreak and not next_moving ) or \
             ( right_before_chainbreak and not prev_moving ) ):
            terminal_res[ name ].append( m )
        if ( ( prev_moving and not next_moving and not right_before_chainbreak ) or \
             ( next_moving and not prev_moving and not right_after_chainbreak ) ):
            extra_min_res[ name ].append( m )

if len (args.extra_flags) > 0:
    if exists( args.extra_flags ):
        extra_flags_benchmark = open( args.extra_flags ).readlines()
    else:
        extra_flags_benchmark = None
        print 'Did not find ', args.extra_flags, ' so not using any extra flags for the benchmark'
        assert ( args.extra_flags == default_extra_flags_benchmark )


for name in names:
    
    dirname = name
    if not exists( dirname ): system( 'mkdir '+dirname )
  
    fid = open( '%s/README_SWM' % name, 'w' )
    fid.write( 'stepwise @flags -out:file:silent swm_rebuild.out\n' )
    fid.close()

    fid = open( '%s/flags' % name, 'w' )

    for infile in [ fasta[name] ] + helix_files[ name ] + input_pdbs[ name ]:  system( 'cp %s %s/ ' % ( infile, name ) )

    start_files = helix_files[ name ] + input_pdbs[ name ]
    if len( start_files ) > 0 :
        fid.write( '-s' )
        for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
        fid.write( '\n' )
    fid.write( '-fasta %s.fasta\n' % name )
    if len( terminal_res[ name ] ) > 0:
        fid.write( '-terminal_res %s  \n' % make_tag_with_conventional_numbering( terminal_res[ name ], resnums[ name ], chains[ name ] ) )
    if len( extra_min_res[ name ] ) > 0 and not args.legacy: ### Turn extra_min_res off for SWM when comparing to SWA
        fid.write( '-extra_min_res %s \n' % make_tag_with_conventional_numbering( extra_min_res[ name ], resnums[ name ], chains[ name ] ) )
    if ( len( input_pdbs[ name ] ) == 0 ):
        fid.write( '-superimpose_over_all\n' ) # RMSD over everything -- better test since helices are usually native
    fid.write( '-cycles 200\n' )
    fid.write( '-nstruct 20\n' )
    fid.write( '-intermolecular_frequency 0.0\n' )
    #fid.write( '-save_times\n' )

    if len( native[ name ] ) > 0:
        system( 'cp %s %s/' % (working_native[name],name) )
        fid.write( '-native %s\n' % basename( working_native[name] ) )
    # case-specific extra flags
    if len( extra_flags[name] ) > 0 : fid.write( '%s\n' % extra_flags[name] )
    # extra flags for whole benchmark
    weights_file = ''
    if extra_flags_benchmark:
        for flag in extra_flags_benchmark:
            if ( flag.find('#') == 0 ): continue
            if ( flag.find( '-score:weights' ) == 0 ): weights_file = string.split( flag )[1]
            fid.write( flag )
        if len( weights_file ) > 0:
            assert( exists( weights_file ) )
            system( 'cp ' + weights_file + ' ' + name )

    fid.close()

    print
    print 'Setting up submission files for: ', name
    CWD = getcwd()
    chdir( name )
    system( 'rosetta_submit.py README_SWM SWM 10 %d' % args.nhours )
    chdir( CWD )

    fid_qsub.write( 'cd %s; source qsubMINI; cd %s\n' % ( name, CWD ) )


    if args.legacy: # SETUP for StepWise Assembly
        
        dirname = name+'/'+'SWA'
        if not exists( dirname ): system( 'mkdir '+dirname )
  
        fid = open( '%s/README_SWA' % dirname, 'w' )
        fid.write( '~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/setup_SWA_RNA_dag_job_files.py' )
        fid.write( ' -single_stranded_loop_mode True' )
        
        if input_res[ name ] != '':
            working_res_ranges = working_res[ name ].split(',')
            input_res_ranges = input_res[ name ].split(',')
            working_res_str = make_tag_from_list_of_int_ranges( working_res_ranges )
            input_res_str = make_tag_from_list_of_int_ranges( input_res_ranges )
            print 'Working_res: '+working_res_str
            print '  Input_res: '+input_res_str
            offset = 1 - int(working_res_str.split()[0])
            sample_res = string.join([ str(int(x)+offset) for x in working_res_str.split() if x not in input_res_str ],' ')     
        else:
            print "ERROR: Must define Input_res for SWA setup."
            exit(0)
        
        fid.write( ' -sample_res %s' % sample_res )

        for infile in [ fasta[name] ] + helix_files[ name ] + input_pdbs[ name ]:  system( 'cp %s %s/ ' % ( infile, dirname ) )

        start_files = helix_files[ name ] + input_pdbs[ name ]
        if len( start_files ) > 0 :
            fid.write( ' -s' )
            for infile in start_files:  fid.write( ' %s' % (basename(infile) ) )
        fid.write( ' -fasta %s.fasta' % name )
        if len( native[ name ] ) > 0:
            system( 'cp %s %s/' % (working_native[name],dirname) )
            fid.write( ' -native_pdb %s' % basename( working_native[name] ) )
        
        # case-specific extra flags
        if len( extra_flags[name] ) > 0 : fid.write( ' %s' % extra_flags[name] )
        # extra flags for whole benchmark
        weights_file = ''
        if extra_flags_benchmark:
            for flag in extra_flags_benchmark:
                if ( flag.find('#') == 0 ): continue
                if ( flag.find( '-score:weights' ) == 0 ):
                    flag.replace( '-score:weights', '-force_field_file' ) 
                    weights_file = string.split( flag )[1]
                if ( flag.find( '-score:rna_torsion_potential' ) == 0 ):
                    flag.replace( '-score:rna_torsion_potential', '-rna_torsion_potential_folder' )
                fid.write( flag )
        
            #if len( weights_file ) > 0:
            #   assert( exists( weights_file ) )
            #   system( 'cp ' + weights_file + ' ' + name )

        fid.close()

        print
        print 'Setting up submission files for: ', name
        CWD = getcwd()
        fid_submit = open( dirname+'/SUBMIT_SWA', 'w' )
        fid_submit.write( '~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/dagman/submit_DAG_job.py' )
        fid_submit.write( ' -master_wall_time %d'% 72 ) #args.nhours )
        fid_submit.write( ' -master_memory_reserve 2048' )
        fid_submit.write( ' -num_slave_nodes 50' )
        fid_submit.write( ' -dagman_file rna_build.dag' )
        fid_submit.close()

        fid_qsub.write( 'cd %s; source ./README_SWA && source ./SUBMIT_SWA; cd %s\n' % ( dirname, CWD ) )

fid_qsub.close()

