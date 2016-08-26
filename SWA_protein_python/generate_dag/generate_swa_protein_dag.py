#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser,expandvars,abspath
from sys import exit, argv
from glob import glob
import string
from time import sleep
from parse_options import parse_options
from get_sequence import get_sequence
from make_tag import make_tag, make_tag_with_dashes

###################
# Clusterer fix
###################
FIX_CALC_RMS_TAG = 1

fasta_file = parse_options( argv, "fasta", "1shf.fasta" )
assert( exists( fasta_file ) )
sequence_lines = open( fasta_file  ).readlines()[1:]
sequence = string.join(  map( lambda x : x[:-1], sequence_lines) ,  '' )
NRES = len( sequence )
NRES = parse_options( argv, "nres_to_model", NRES )
EXE = parse_options( argv, "exe", "" )
DB = parse_options( argv, "database", "" )
MIN_RES = parse_options( argv, "min_res", 1 )
MAX_RES = parse_options( argv, "max_res", NRES )
ZIGZAG = parse_options( argv, "zigzag", 0 )
N_SAMPLE = parse_options( argv, "n_sample", 18 )
FINAL_NUMBER = parse_options( argv, "final_number", 100 )
SCORE_WEIGHTS = parse_options( argv, "weights", "score12_no_hb_env_dep.wts" )
PACK_WEIGHTS = parse_options( argv, "pack_weights", "pack_no_hb_env_dep.wts" )
NSTRUCT = parse_options( argv, "nstruct", 400 )
RMSD_SCREEN = parse_options( argv, "rmsd_screen", -1.0 )
CLUSTER_RADIUS = parse_options( argv, "cluster_radius", 0.25 )
CLUSTER_RADIUS_SAMPLE = parse_options( argv, "cluster_radius_sample", 0.1 )
AUTO_TUNE = parse_options( argv, "auto_tune", 0 )
filter_native_big_bins = parse_options( argv, "filter_native_big_bins", 0 )
move_all_jumps = parse_options( argv, "move_all_jumps", 0 )
move_jumps_between_chains = parse_options( argv, "move_jumps_between_chains", 0 )
score_diff_cut = parse_options( argv, "score_diff_cut", 10.0 )
max_res_to_add_denovo = parse_options( argv, "denovo", 0 )
USE_MINI_TEMP = parse_options( argv, "use_mini_TEMP", 0 )
no_rm_files= parse_options( argv, "no_rm_files", 0 )
secstruct= parse_options( argv, "secstruct", "" )

native_pdb = parse_options( argv, "native", "" )
template_pdbs = parse_options( argv, "s", [""] )
template_input_res = parse_options( argv, "input_res", [-1] )
cst_file = parse_options( argv, "cst_file", "" )
disulfide_file = parse_options( argv, "disulfide_file", "" )
pathway_file = parse_options( argv, "pathway_file", "" )
cluster_by_all_atom_rmsd = parse_options( argv, "cluster_by_all_atom_rmsd", 0 )
add_peptide_plane = parse_options( argv, "add_peptide_plane", 0 ) #Now defunct!
no_peptide_plane = parse_options( argv, "no_peptide_plane", 0 )
BUILD_BOTH_TERMINI = parse_options( argv, "build_both_termini", 0 )
MAX_ADDED_SEGMENT = parse_options( argv, "max_added_segment", 20 )
min_length = parse_options( argv, "min_length", 2 )
max_length = parse_options( argv, "max_length", 0 )
superimpose_res = parse_options( argv, "superimpose_res", [ -1 ] )
virtual_res = parse_options( argv, "virtual_res", [ -1 ] )
skip_res = parse_options( argv, "skip_res", [ -1 ] )
jump_res = parse_options( argv, "jump_res", [ -1 ] )
align_pdb = parse_options( argv, "align_pdb", "" )
template_mapping_files = parse_options( argv, "mapping", [""] )
frag_files = parse_options( argv, "frag_files", [""] )
frag_lengths = parse_options( argv, "frag_lengths", [ -1 ] )
MAX_FRAGMENT_OVERLAP = parse_options( argv, "max_fragment_overlap", 2 )
swa_frag_lengths = parse_options( argv, "swa_frag_lengths", [-1] )
swa_silent_file_dir = parse_options( argv, "swa_silent_file_dir", "" )
fixed_res = parse_options( argv, "fixed_res", [-1] )
no_fixed_res = parse_options( argv, "no_fixed_res", 0 )
calc_rms_res = parse_options( argv, "calc_rms_res", [-1] )
start_res = parse_options( argv, "start_res", [-1] )
loop_start_pdb = parse_options( argv, "loop_start_pdb", "" )
start_pdb = parse_options( argv, "start_pdb", loop_start_pdb )
if len( start_pdb ) > 0:
    start_pdbs = [ start_pdb ]
elif len( loop_start_pdb ) > 0:
    start_pdbs = [ loop_start_pdb ]
else: start_pdbs = parse_options( argv, "start_pdbs", [""] )

loop_start_full_outfile = parse_options( argv, "loop_start_full_outfile", "" )
endpoints = parse_options( argv, "endpoints",[-1] )
loop_res = parse_options( argv, "loop_res", [-1] )
loop_force_Nsquared = parse_options( argv, "loop_force_Nsquared", 0 )
override = parse_options( argv, "override", 0 )
new_loop_close = parse_options( argv, "new_loop_close", 0 )
cutpoints_open = parse_options( argv, "cutpoint_open", [ -1 ] )
cutpoints_closed  = parse_options( argv, "cutpoint_closed", [ -1 ] )
only_go_midway = parse_options( argv, "only_go_midway", 0 )
centroid = parse_options( argv, "centroid", 0 )
DO_CCD = not parse_options( argv,'no_ccd',0 )
DO_KIC = parse_options( argv, 'do_kic', 0 )
disable_sampling_of_loop_takeoff = parse_options( argv, 'disable_sampling_of_loop_takeoff', 0 )

if ( len( argv ) > 1 ): # Should remain with just the first element, the name of this script.
    print " Unrecognized flags?"
    print "   ",string.join(argv[1:] )
    exit( 0 )

DENOVO = ( max_res_to_add_denovo > 0 )
TEMPLATE = len( template_pdbs ) > 0
FRAGMENT_LIBRARY = len( frag_files ) > 0
SWA_FRAGS = len( swa_frag_lengths ) > 0

LOOP = 0
if len( loop_start_pdb ) > 0:
    LOOP = 1
    start_pdb = loop_start_pdb
if len( loop_start_full_outfile ) > 0:
    assert( len( loop_start_pdb) == 0 )
    assert( loop_start_full_outfile == 'region_0_0_sample.cluster.out' )
    LOOP = 1


for template_pdb in template_pdbs: assert( exists( template_pdb ) )
for frag_file in frag_files: assert( exists( frag_file ) )
if add_peptide_plane:
    print " -add_peptide_plane defunct -- its on by default! "
    print " If you want to disable peptide_plane, use no_peptide_plane."
    exit( 0 )
add_peptide_plane = not no_peptide_plane
if FRAGMENT_LIBRARY: assert( add_peptide_plane )   # Should peptide plane also be required for template runs?

###############################################################
# Where's the executable?
###############################################################

#MINI = "mini"
MINI = "rosetta_TRUNK/rosetta_source/"
if USE_MINI_TEMP: MINI = "mini_TEMP"

rosetta_folder = ''
if len( EXE ) == 0:
    rosetta_folder = expandvars("$ROSETTA")
    if rosetta_folder == "$ROSETTA" :
        error_exit("USERs need to set environmental variable $ROSETTA and pointed it to the Rosetta folder!")
    exe_folder = rosetta_folder + "/rosetta_source/bin/" #Default Rosetta folder structure
    if not exists(exe_folder) : #Otherwise, assume the input folder name is bin path
        exe_folder = rosetta_folder
    name_extensions = [".linuxgccrelease", ".linuxclangrelease", ".macosgccrelease", ".macosclangrelease",
                       "failed_to_find_Rosetta_path"] #this makes a better error message if pathing fails
    exe_file = "swa_protein_main"
    exe_path = ""
    for name in name_extensions :
        EXE = exe_folder + exe_file + name
        if exists(EXE) :
            break

assert( exists( EXE ) )

if len( DB ) == 0 and len( rosetta_folder ) > 0:
    DB = rosetta_folder + "/rosetta_database/"
assert( exists( DB ) )

PYDIR = abspath( dirname(argv[0]) + "/../run_dag_on_cluster/" )
print PYDIR
assert( exists( PYDIR ) )

PRE_PROCESS_SETUP_SCRIPT = PYDIR+"/stepwise_pre_process_setup_dirs.py"
assert( exists( PRE_PROCESS_SETUP_SCRIPT ) )

POST_PROCESS_FILTER_SCRIPT = PYDIR+"/stepwise_post_process_combine_and_filter_outfiles.py"
assert( exists( POST_PROCESS_FILTER_SCRIPT ) )

POST_PROCESS_CLUSTER_SCRIPT = PYDIR+"/stepwise_post_process_cluster.py"
assert( exists( POST_PROCESS_CLUSTER_SCRIPT ) )

assert( exists( SCORE_WEIGHTS ) or exists( DB + "/scoring/weights/"+SCORE_WEIGHTS) )
assert( exists( PACK_WEIGHTS ) or exists( DB + "/scoring/weights/"+PACK_WEIGHTS) )

fid_dag = open( "protein_build.dag", 'w' )
fid_dag.write("DOT dag.dot\n")

if no_rm_files:
    POST_PROCESS_CLUSTER_SCRIPT += ' -no_rm_files 1'
    POST_PROCESS_FILTER_SCRIPT += ' -no_rm_files 1'

if not exists( 'CONDOR/' ):
    system( 'mkdir -p CONDOR' )
    #system( 'chmod 777 -R CONDOR' )

#########################################################
# list of jobs...
all_job_tags = []
real_compute_job_tags = []
jobs_done = []

# Special case.
if len( loop_start_full_outfile ) > 0 :
    all_job_tags.append( 'REGION_0_0' )
    jobs_done.append( 'REGION_0_0' )

#########################################################
# Some useful functions (move somewhere else?)
#########################################################
def make_condor_submit_file( condor_submit_file, arguments, queue_number, universe="vanilla" ):

    fid = open( condor_submit_file, 'w' )
    fid.write('+TGProject = TG-MCB090153\n')
    fid.write('universe = %s\n' % universe)
    fid.write('executable = %s\n' % EXE )

    fid.write('arguments = %s\n' % arguments)

    sub_job_tag = basename( condor_submit_file ).replace('.condor','')
    job_dir = dirname( condor_submit_file )
    sub_job_dir = job_dir + '/' + sub_job_tag

    assert( exists( job_dir ) )
    if not exists( sub_job_dir ):
        system( 'mkdir -p '+sub_job_dir )
        ##system( 'chmod 777 -R '+sub_job_dir )

    fid.write('output = %s/$(Process).out\n' % sub_job_dir )
    fid.write('log = %s/%s.log\n' % ( job_dir,sub_job_tag) )
    fid.write('error = %s/$(Process).err\n' % sub_job_dir )
    fid.write('notification = never\n')
    fid.write('Queue %d\n' % queue_number )
    fid.close()

def setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag,\
                                         fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files):
    newdir = overall_job_tag+'/'+sub_job_tag
    if len( decoy_tag ) > 0:
        newdir += '_' + decoy_tag
    outfile = newdir + '/' + overall_job_tag.lower() + '_sample.out'

    args2 += ' -out:file:silent %s ' % outfile

    condor_file_dir =  "CONDOR/%s/" % overall_job_tag
    if not exists( condor_file_dir ):
        system( 'mkdir -p '+condor_file_dir )
        ##system( 'chmod 777 -R '+condor_file_dir )

    condor_submit_file = '%s/%s.condor' %  (condor_file_dir,sub_job_tag)

    job_tag = overall_job_tag+"_"+sub_job_tag
    fid_dag.write('\nJOB %s %s\n' % (job_tag, condor_submit_file) )

    if not exists( condor_submit_file ):
        make_condor_submit_file( condor_submit_file, args2, 1 )

    if (len( prev_job_tags ) > 0):
        for prev_job_tag in prev_job_tags:
            assert( prev_job_tag in all_job_tags )
            #Note previous job may have been accomplished in a prior run -- not in the current DAG.
            if (prev_job_tag not in jobs_done):
                fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, job_tag) )

    # The pre process script finds out how many jobs there actually are...
    if len( decoy_tag ) > 0:
        assert( len( prev_job_tags ) > 0 )
        fid_dag.write('SCRIPT PRE %s   %s %s %s %s %s\n' % (job_tag, PRE_PROCESS_SETUP_SCRIPT,overall_job_tag,prev_job_tags[0],condor_submit_file,sub_job_tag) )
    else:
        if not exists( newdir ):
            system( 'mkdir -p '+newdir )
            ##system( 'chmod 777 -R '+newdir )

    fid_dag.write('SCRIPT POST %s %s %s/%s\n' % (job_tag, POST_PROCESS_FILTER_SCRIPT,overall_job_tag,sub_job_tag ) )

    # In python, lists are passed by reference... these should get updated for the outside world.
    job_tags.append( job_tag )
    real_compute_job_tags.append( job_tag )
    combine_files.append( '%s/%s_sample.low4000.out' % ( overall_job_tag, sub_job_tag.lower() ) )


#####################################################
# Pathway setup
#####################################################
follow_path = 0
if len(pathway_file) > 0:
    follow_path = 1
    lines = open( pathway_file ).readlines()
    pathway_regions = []

    for line in lines:
        # don't yet know how to handle "merges" (i.e., inter-domain docking)
        if len( line ) > 5  and line[:5] == 'MERGE': break
        assert( line[:4] == 'PATH' )

        cols = map( lambda x:int(x), string.split( line )[1:] )
        i = cols[0]
        j = cols[0]
        cols = cols[1:]
        for m in cols:
            if ( m == i-1 ):
                i = m
            else:
                assert( m == j+1 )
                j = m
            pathway_regions.append( [i, j] )

#####################################################
# Template mapping files
#####################################################
if len( template_pdbs ) > 0:
    template_mappings = {}
    if len( template_mapping_files ) == 0:
        template_mapping = {}
        if len( template_input_res ) > 0:
            count = 0
            check_sub_sequence = ''
            for m in template_input_res:
                count += 1
                template_mapping[ m ] =  count
                check_sub_sequence += sequence[ m-1 ]
            template_sequence = get_sequence( template_pdbs[ 0 ] )
            #print check_sub_sequence
            #print template_sequence
            assert( check_sub_sequence == template_sequence )
        else:
            # Asssume that sequences correspond perfectly.
            for m in range( 1, NRES+1 ): template_mapping[ m ] = m
        for n in range( len( template_pdbs) ):
            template_mappings[ template_pdbs[n] ] =  template_mapping
    else:
        assert( len( template_mapping_files ) == len( template_pdbs ) )
        for n in range( len( template_pdbs ) ):
            lines = open( template_mapping_files[n] ).readlines()
            # First line better be our target sequence
            mapping_seq1 = lines[0][:-1]
            #print mapping_seq1.replace('-','')
            #print sequence
            assert(  mapping_seq1.replace('-','') == sequence )

            # Second line better be our template sequence
            template_sequence = popen( "python "+PYDIR+"/pdb2fasta.py "+template_pdbs[n] ).readlines()[1][:-1]
            mapping_seq2 = lines[1][:-1]
            #print mapping_seq2.replace('-','')
            #print sequence
            assert( mapping_seq2.replace('-','') == template_sequence )
            assert( len( mapping_seq1 ) == len( mapping_seq2 ) )

            count1 = 0
            count2 = 0
            template_mapping = {}
            for i in range( len( mapping_seq1 ) ):
                seq1_OK = ( not mapping_seq1[i]  == '-' )
                seq2_OK = ( not mapping_seq2[i]  == '-' )
                if seq1_OK: count1 += 1
                if seq2_OK: count2 += 1
                if seq1_OK and seq2_OK: template_mapping[ count1 ] = count2
            template_mappings[ template_pdbs[n] ] = template_mapping

def template_continuous( i, j, template_mapping ):
    for k in range( i, j+1 ):
        if not ( k in template_mapping.keys() ): return 0
        if ( k > i ) and ( not template_mapping[k] == template_mapping[k-1]+1 ): return 0
    return 1

##########################
# Fragments
##########################

if len( frag_lengths ) == 0:
    for frag_file in frag_files:
        pos = frag_file.find("_05.200_v1_3" )
        assert( pos > 0 )
        frag_length = int( frag_file[pos-2 : pos] )
        frag_lengths.append( frag_length )

assert( len( frag_lengths ) == len( frag_files ) )


def check_frag_overlap( i, j, i_prev, j_prev, frag_length ):

    input_res1 =  wrap_range( i_prev, j_prev+1 )

    if ( i == i_prev ):
        assert( not j == j_prev )
        num_extra_residues = j - j_prev
        startpos = j - frag_length + 1
    elif ( j == j_prev ):
        num_extra_residues = i_prev - i
        startpos = i
    else:
        # Some craziness
        return (0,0,0)

    num_overlap_residues = frag_length - num_extra_residues

    endpos = startpos + frag_length - 1

    return ( startpos, endpos, num_overlap_residues )

####################################################################
# input silent files from some other directory (e.g., SWA frags)
####################################################################
if SWA_FRAGS:
    glob_files = glob( swa_silent_file_dir+"/region*sample.cluster.out" )
    silent_files_in = map( lambda x: basename(x),  glob_files )

#####################################
# Setup for loop building jobs
#####################################
def wrap_range( i, j, total_residues = NRES ):
    # Kind of like range(i,j).
    # But if i > j-1, gives residues not in j ... i-1.
    # Useful for loop stuff
    new_range = []

    if ( i < j ):
        new_range = range(i,j)
    else:
        new_range = range(1,j)
        for m in range( i, total_residues+1 ): new_range.append( m )

    return new_range

def wrap_range_and_add_satellite( i, j, total_residues = NRES, start_res_for_input_pdb = []):
    new_range = wrap_range( i, j, NRES )


    for n in range( len(start_res_for_input_pdb) ):

        start_res = start_res_for_input_pdb[ n ]

        found_a_start_res = 0
        for m in start_res:
            if m in new_range:
                found_a_start_res = 1
                break

        if found_a_start_res:
            for m in start_res:
                if m not in new_range:
                    new_range.append( m )

        new_range.sort()

    return new_range

def reorder_to_be_contiguous( res_list_input, total_residues = NRES ):
    in_order = 1

    res_list = []
    for m in res_list_input:
        if (m <= NRES): res_list.append( m )

    for n in range( len( res_list )-1 ):
        if not ( res_list[n+1]-1 == res_list[n] ):
            in_order = 0
            break
    if in_order:
        return res_list
    else:
        res_list_new = []
        for m in range( n+1, len(res_list) ):  res_list_new.append( res_list[m] )
        for m in range( n+1 ):                 res_list_new.append( res_list[m] )
        for n in range( len(res_list_new)-1 ):
            assert( res_list_new[n+1]-1 == res_list_new[n] or ( res_list_new[n] == total_residues and res_list_new[n+1] == 1 ) )
        return res_list_new

min_loop_gap = 2
if ( new_loop_close ): min_loop_gap = 4

START_FROM_PDB = len( start_pdbs ) > 0  or len( loop_start_full_outfile ) > 0
if len( loop_start_full_outfile ) > 0: assert( exists( loop_start_full_outfile ) )
if LOOP: assert( START_FROM_PDB )

start_res_for_input_pdb = []
cutpoint_open_in_loop = 0

all_loop_res = []
LOOP_AT_TERMINUS = 0
if START_FROM_PDB:

    # -loop_res for one contiguous loop was specified, not -start_res
    # This was the main mode used for a lot of CASP9 rebuilds.
    # However, I think it would be better to define everything in terms
    # of -start_res (i.e., the stuff that isn't loop), since this
    # will allow us to deal with topologies in which there are multiple interacting loops.

    if len( loop_res ) > 0:
        assert( len(start_res) == 0 )
        loop_res.sort()
        assert( len( loop_res ) >= min_loop_gap)
        for m in range( len(loop_res)-1 ): assert( loop_res[m]+1 == loop_res[m+1] )

        loop_start = loop_res[ 0 ]
        loop_end = loop_res[ -1 ]

        if loop_end == NRES or loop_start == 1: LOOP_AT_TERMINUS = 1

        start_res = wrap_range( loop_end+1, loop_start, NRES )
        #start_res_for_input_pdb.append( start_res )

        if len( calc_rms_res ) == 0:  calc_rms_res = range( loop_start, loop_end+1 )
        if len( superimpose_res ) == 0: superimpose_res = start_res

        if len( endpoints ) > 0 and ( loop_start-1 not in endpoints ):
            obligate_endpoints = [ loop_start-1, loop_start, loop_start+1, loop_start+2, \
                                   loop_end-2, loop_end-1, loop_end, loop_end+1 ]
            for m in obligate_endpoints:
                if m not in endpoints:
                    print "Adding ", m, " to endpoints"
                    endpoints.append( m )

        all_loop_res.append( loop_res )

    else:
        if ( start_res ) == 0:
            # This was also used in CASP9, but probably should be deprecated.
            # It makes much better sense to have the user input the start_res manually,
            # and then make sure the sequence looks OK in here. Consistency checks!!
            assert( len ( start_pdbs ) == 1)
            assert( exists( start_pdb ) )
            start_sequence = get_sequence( start_pdb )
            nres_start = len( start_sequence )
            found_match = 0
            for i in range( NRES - nres_start + 1):
                if (sequence[ i:(i+nres_start) ] == start_sequence):
                    found_match = 1
                    break
            assert( found_match )
            start_res = range( i+1, (i+nres_start+1) )

            #if FIX_CALC_RMS_TAG and len( calc_rms_res ) == 0:
            if ( i > 1 ):
                calc_rms_res = range( 1, i+1 )
            else:
                assert( i+nres_start+1 < NRES )
                calc_rms_res = range( i+nres_start+1, NRES+1 )

            print "Figured out that starting pdb %s has residues %d-%d" % ( start_pdb, i+1, i+nres_start )

    ##############################################################
    # Consistency check and split of user-defined start_res
    #  across multiple input pdbs...
    assumed_start_sequence = ''
    #print start_res
    for m in start_res:
        assumed_start_sequence += sequence[m-1]

    full_actual_start_sequence = ''
    count = 0
    for start_pdb in start_pdbs:
        actual_start_sequence = get_sequence( start_pdb )
        start_res_for_this_pdb = []
        for m in actual_start_sequence:
            full_actual_start_sequence += m
            start_res_for_this_pdb.append( start_res[count] )
            count += 1
            if ( 1 in start_res_for_this_pdb and NRES in start_res_for_this_pdb ):
                jump_res.append( 1 )
                jump_res.append( NRES )

        #print  start_res_for_this_pdb, reorder_to_be_contiguous( start_res_for_this_pdb )
        #start_res_for_input_pdb.append( reorder_to_be_contiguous( start_res_for_this_pdb ) )
        start_res_for_input_pdb.append( start_res_for_this_pdb )

    #print assumed_start_sequence, len( assumed_start_sequence)
    #print full_actual_start_sequence, len( full_actual_start_sequence )
    assert( assumed_start_sequence == full_actual_start_sequence )
    #################################

    if len( endpoints ) > 0 and not LOOP:
        obligate_endpoints = [start_res[0], start_res[-1]]
        for m in obligate_endpoints:
            if m not in endpoints:
                print "Adding ", m, " to endpoints"
                endpoints.append( m )

    if len( start_res_for_input_pdb ) > 1:
        print "User is asking for multiple loops. Forcing O(N^2) run"
        loop_force_Nsquared = 1

    # OK, loop residue segments.
    if len( all_loop_res ) == 0:
        this_loop_res = []
        for m in range( 1, NRES+1 ):
            in_start_res = 0
            for start_res in start_res_for_input_pdb:
                if ( m in start_res ):
                    in_start_res = 1
                    break
            if (in_start_res):
                if ( len(this_loop_res) > 0 ): all_loop_res.append( this_loop_res )
                this_loop_res = []
            else:
                this_loop_res.append( m )
    #print all_loop_res

    min_start_pdb_length = NRES
    for start_res in start_res_for_input_pdb:
        start_pdb_length = len( reorder_to_be_contiguous( start_res ) )
        if ( start_pdb_length < min_start_pdb_length ): min_start_pdb_length = start_pdb_length

    if len( fixed_res ) == 0 and not no_fixed_res:
        fixed_res = []
        for start_res in start_res_for_input_pdb:
            for m in start_res: fixed_res.append( m )

    if len( calc_rms_res ) == 0:
        for i in range(1, NRES+1):
            if i not in fixed_res: calc_rms_res.append( i )

    if len( superimpose_res ) == 0: superimpose_res = fixed_res


for k in virtual_res:
    if k not in skip_res:  skip_res.append( k )
skip_res.sort()

for k in skip_res:
    if (k not in fixed_res):  fixed_res.append( k )
fixed_res.sort()

if centroid:
    PACK_WEIGHTS =  'score3_with_cst.wts'
    SCORE_WEIGHTS = 'score3_min.wts'

##########################
# BASIC COMMAND
##########################

args = ' -database %s  -rebuild -out:file:silent_struct_type binary  -fasta %s -n_sample %d -nstruct %d -cluster:radius %8.3f' % ( DB, fasta_file, N_SAMPLE, NSTRUCT, CLUSTER_RADIUS_SAMPLE )

args += ' -extrachi_cutoff 0 -ex1 -ex2' # These may be redundant actually.

args += ' -score:weights %s -pack_weights %s' % (SCORE_WEIGHTS, PACK_WEIGHTS )

 # prevents the occasional automatic disulfide check from ruining things.
 # only disulfides that show up will be the ones specified by -disulfide_file
args += ' -in:detect_disulf false'

if ( RMSD_SCREEN > 0.0 ):
    args += ' -rmsd_screen %8.3f' % RMSD_SCREEN

if add_peptide_plane: args += ' -add_peptide_plane'
if filter_native_big_bins:  args+= ' -filter_native_big_bins' # this is defunct now, I think
if move_all_jumps:  args+= ' -move_all_jumps'
if move_jumps_between_chains:  args+= ' -move_jumps_between_chains'

if len( cst_file ) > 0:
    assert( exists( cst_file ) )
    args += ' -cst_file %s' % cst_file
if len( disulfide_file ) > 0:
    assert( exists( disulfide_file ) )
    args += ' -disulfide_file %s' % disulfide_file
if len( align_pdb ) > 0:
    assert( exists( align_pdb ) )
    args += ' -align_pdb %s' % align_pdb
if len( native_pdb ) > 0:
    assert( exists( native_pdb ) )
    args += ' -native %s' % native_pdb
if len( superimpose_res ) > 0:
    args += ' -superimpose_res '
    args += make_tag_with_dashes( superimpose_res )
if len( virtual_res ) > 0:
    args += ' -virtual_res '
    args += make_tag_with_dashes( virtual_res )
if len( fixed_res ) > 0:
    args += ' -fixed_res '
    args += make_tag_with_dashes( fixed_res )
if len( calc_rms_res ) > 0:
    args += ' -calc_rms_res '
    args += make_tag_with_dashes( calc_rms_res )
if len( jump_res ) > 0:
    args += ' -jump_res '
    args += make_tag( jump_res )
if len( cutpoints_open ) > 0:
    args += ' -cutpoint_open '
    args += make_tag_with_dashes( cutpoints_open )
if len( cutpoints_closed ) > 0:
    args += ' -cutpoint_closed '
    args += make_tag_with_dashes( cutpoints_closed )
if centroid:
    args += ' -centroid '
    #args += ' -centroid -skip_minimize '
if len( secstruct ) > 0:
    args += ' -secstruct '+secstruct
if disable_sampling_of_loop_takeoff:
    args += ' -disable_sampling_of_loop_takeoff'

def add_cutpoint_closed( args, cutpoint ):
    if args.find( '-cutpoint_closed' ) > -1:
        args = args.replace( '-cutpoint_closed', '-cutpoint_closed %d' % cutpoint )
    else:
        args += ' -cutpoint_closed %d' % cutpoint
    return args

args += ' -mute all' # trying to cut down on disk space!

args_START = args

if AUTO_TUNE:
    cluster_tag = ' -auto_tune '
else:
    cluster_tag = ' -cluster:radius %s ' % CLUSTER_RADIUS


if FIX_CALC_RMS_TAG and len( calc_rms_res ) > 0:
    cluster_tag += ' -calc_rms_res'
    for k in calc_rms_res: cluster_tag += ' %d' % k

cluster_by_all_atom_rmsd_tag = ''
if cluster_by_all_atom_rmsd: cluster_by_all_atom_rmsd_tag = ' -cluster_by_all_atom_rmsd '

def check_in_another_loop( i, loop_start, all_loop_res ):
    for loop_res in all_loop_res:
        if loop_res[0] == loop_start: continue
        if i in loop_res: return True
    return False


################################
# MAIN LOOP
################################
# Loop over fragment lengths.

if (max_length == 0): max_length = NRES
for L in range( min_length, max_length + 1 ):
    chunk_length = L;
    #num_chunks = ( len( sequence) - chunk_length) + 1

    if START_FROM_PDB and (L < min_start_pdb_length): continue

    for k in range( 1, NRES + 1 ) :
        i = k
        j = i + chunk_length - 1
        if ( j > NRES ): j -= NRES
        res_to_be_modeled = wrap_range(i,j+1,NRES)

        if ( ( not START_FROM_PDB ) and ( i < MIN_RES or j > MAX_RES or i > j ) ): continue
        if ( i in skip_res or j in skip_res ): continue

        loop_close = 0

        if START_FROM_PDB:

            if LOOP_AT_TERMINUS and i>j: continue

            #######################################
            # Do we contain one of the starting PDBs?
            contains_at_least_one_start_region = 0
            does_not_contain_whole_start_region = 0
            for start_res in start_res_for_input_pdb:
                contains_start_region = 0
                for m in reorder_to_be_contiguous(start_res):
                    if m in res_to_be_modeled:
                        contains_start_region = 1
                        break

                if contains_start_region:
                    contains_at_least_one_start_region = 1
                    for m in reorder_to_be_contiguous(start_res):
                        if m not in res_to_be_modeled:
                            does_not_contain_whole_start_region = 1
                            break
                    if does_not_contain_whole_start_region: break

            if ( not contains_at_least_one_start_region ): continue
            if ( does_not_contain_whole_start_region ): continue

            #######################
            # Setup for loops...
            #######################
            # which loop are we in?
            found_specific_loop = 0
            for loop_res in all_loop_res:

                loop_res_expand = [ loop_res[0]-1 ]
                for m in loop_res: loop_res_expand.append( m )
                loop_res_expand.append( loop_res[-1]+1 )

                if ( (j in loop_res_expand or (loop_res[-1] == NRES and j == 1) ) and
                     (i in loop_res_expand or (loop_res[ 0] == 1 and i == NRES) ) ):
                    found_specific_loop = 1
                    loop_start = loop_res[0]
                    loop_end = loop_res[-1]
                    break

            ###########################################
            loop_close = (i == j+1)
            if loop_close and not found_specific_loop: continue
            if ( i == 1 and j == NRES and not found_specific_loop and not LOOP_AT_TERMINUS): continue

            if found_specific_loop:
                # Don't build across open cutpoints
                found_cutpoint_open = 0
                for m in range( loop_start-1, j):
                    if m in cutpoints_open:
                        found_cutpoint_open = 1
                        break
                for m in range( i, loop_end+1):
                    if m in cutpoints_open:
                        found_cutpoint_open = 1
                        break

                if found_cutpoint_open: continue


                # Even if we don't build across cutpoints, this loop may have an open cutpoint later...
                # That determine whether or not we bother building all the way to the end
                cutpoint_open_in_loop = 0
                for m in range( loop_start-1, loop_end+1 ):
                    if m in cutpoints_open: cutpoint_open_in_loop = 1
                if ( i == 1 or  j == NRES ): cutpoint_open_in_loop = 1
                if ( cutpoint_open_in_loop ): loop_close = 0


                # How far forward or backward should we build loop?
                # Can't go all the way to the end -- we need to leave a gap for loop closure.
                #
                #  The one special case here is if there is a cutpoint_open at the loop boundary
                #  we won't close the chain -- we'll just build to the end.
                #if ( not cutpoint_open_in_loop  and len( all_loop_res ) == 1 and i < (loop_start+2)   ): continue
                #if ( not cutpoint_open_in_loop  and len( all_loop_res ) == 1 and j > (loop_end  -2)   ): continue

                # boundary cases that do not proceed into loop closure
                #if ( not cutpoint_open_in_loop and i == loop_end+1     and j == loop_end-1   ): continue
                #if ( not cutpoint_open_in_loop and i == loop_start+1   and j == loop_start-1 ): continue

                if not cutpoint_open_in_loop:
                    if loop_close:
                        #if j < loop_start: continue
                        #if i > loop_end-1: continue
                        if j < loop_start-1: continue
                        if i > loop_end+1: continue
                    else:
                        #if (i-j) <= 3: continue
                        if (i-j) <= 1: continue

                # very special case for building termini -- this is now covered above
                #if ( loop_start == 1  and j == 1 ): continue
                #if ( loop_end == NRES and i == NRES ): continue
                #if not( ( loop_close and j < (loop_end-1) )   or ( (i - j) >= 2 ) ): continue

                # To close loop, must have at least a little bit built from either end.
                #if ( loop_close and ( i >= loop_end or j < loop_start ) and not cutpoint_open_in_loop ): continue
                # actually in new scheme, do not need anything prebuilt.

                # Unless special O(N^2) type run (sample little bits of loop in both forward and reverse directions ),
                #  force either i or j to be at boundary.
                if (not loop_force_Nsquared) and  (not loop_close) and not ( i == loop_end+1 or j == loop_start-1): continue

                if ( only_go_midway ):
                    if ( i == loop_end+1    and  j > (( loop_start+loop_end)/2 + 1) ): continue
                    if ( j == loop_start-1  and  i < (( loop_start+loop_end)/2 - 1) ): continue
                    if ( loop_close    and  j > (( loop_start+loop_end)/2 + 1) ): continue
                    if ( loop_close    and  i < (( loop_start+loop_end)/2 - 2) ): continue

                if ( loop_close and disable_sampling_of_loop_takeoff ):
                    if ( i == loop_end + 1): continue
                    if ( j == loop_start - 1): continue


        if len( endpoints ) > 0:
            if ( i not in endpoints ): continue
            if ( j not in endpoints ): continue

        #ZIGZAG!! special case for beta hairpins.
        if ( ZIGZAG and abs( ( i - MIN_RES ) - ( MAX_RES - j ) ) > 1 ) : continue

        if follow_path and ( [i,j] not in pathway_regions ): continue

        overall_job_tag = 'REGION_%d_%d' % (i,j)

        print 'Do region ==> %3d %3d   ' %(i,j) ,

        # This job is maybe already done...
        outfile_cluster = overall_job_tag.lower()+'_sample.cluster.out'
        if exists( outfile_cluster ):
            all_job_tags.append(  overall_job_tag )
            jobs_done.append( overall_job_tag   )
            print 'DONE'
            continue

        termini_tag = ""
        if ( i == 1 ): termini_tag += " -n_terminus"
        if ( j == NRES ): termini_tag += " -c_terminus"
        args = args_START + termini_tag


        ###########################################
        # DO THE JOBS
        ###########################################
        start_regions = []


        for k in range( 1, MAX_ADDED_SEGMENT):
            i_prev = i
            j_prev = j - k
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
            if prev_job_tag in all_job_tags:   start_regions.append( [i_prev, j_prev ] )

        for k in range( 1, MAX_ADDED_SEGMENT):
            i_prev = i + k
            j_prev = j
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
            if prev_job_tag in all_job_tags:   start_regions.append( [i_prev, j_prev ] )

        # Is this deprecated?
        if BUILD_BOTH_TERMINI:
            i_prev = i + 1
            j_prev = j - 1
            prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
            if prev_job_tag in all_job_tags:   start_regions.append( [i_prev, j_prev ] )

        # Another possibility -- if we contain a long segment of virtual residues.
        skip_sub_segment = []
        for m in wrap_range(i,j+1):
            if m in skip_res: skip_sub_segment.append( m )
        if len( skip_sub_segment ) > 0 and ( not loop_close ):
            skip_sub_segment.sort()
            skip_start = skip_sub_segment[0]
            skip_end   = skip_sub_segment[-1]
            if (j > skip_end ) and ( j - skip_end ) <= MAX_FRAGMENT_OVERLAP:
                i_prev = i
                for j_prev in wrap_range( skip_start - MAX_FRAGMENT_OVERLAP, skip_start+1 ):
                    if ( [i_prev, j_prev ] in start_regions ): continue
                    prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
                    if prev_job_tag in all_job_tags: start_regions.append( [i_prev, j_prev ] )
                    print ' -- long frag to bridge over skip segment. start from:  %d-%d' % (i_prev,j_prev)
            if (i < skip_start ) and ( skip_start-i ) <= MAX_FRAGMENT_OVERLAP:
                j_prev = j
                for i_prev in wrap_range( skip_end+1, skip_end+MAX_FRAGMENT_OVERLAP+1):
                    if ( [i_prev, j_prev ] in start_regions ): continue
                    prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
                    if prev_job_tag in all_job_tags: start_regions.append( [i_prev, j_prev ] )
                    print ' -- long frag to bridge over skip segment. start from:  %d-%d' % (i_prev,j_prev)

        # One final possibility. Addition of pre-existing rigid chunk (a "starting pdb" )
        if START_FROM_PDB:
            for start_res_check in start_res_for_input_pdb:
                start_res = reorder_to_be_contiguous( start_res_check )
                if ( start_res[0] == i and start_res[-1] in res_to_be_modeled ):
                    i_prev = start_res[-1] + 1
                    j_prev = j
                    prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
                    if prev_job_tag in all_job_tags  and  [i_prev,j_prev] not in start_regions:   start_regions.append( [i_prev, j_prev ] )
                if ( start_res[0] in res_to_be_modeled and start_res[-1] == j):
                    i_prev = i
                    j_prev = start_res[0] - 1
                    prev_job_tag = 'REGION_%d_%d' % (i_prev,j_prev)
                    if prev_job_tag in all_job_tags  and  [i_prev,j_prev] not in start_regions:   start_regions.append( [i_prev, j_prev ] )

        job_tags = []
        combine_files = []

        ########################################################
        # Several "modes" -- different stuff for the grinder!
        #
        #  currently have: DENOVO, TEMPLATE, FRAGMENT_LIBRARY
        ########################################################
        prev_job_tags = []
        if START_FROM_PDB:

            for n in range( len(start_res_for_input_pdb) ):
                start_res = start_res_for_input_pdb[ n ]


                if ( reorder_to_be_contiguous( res_to_be_modeled ) == reorder_to_be_contiguous( start_res ) ):

                    if ( len( loop_start_full_outfile ) > 0 ):
                        sub_job_tag = 'START_FROM_OUTFILE'

                        decoy_tag = 'S_$(Process)'

                        infile1 =  loop_start_full_outfile
                        args2 = '%s  -silent1 %s -tags1 %s' % (args, infile1, decoy_tag )

                        args2 += " -input_res1 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )

                        args2 += " -slice_res1 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )

                        args2 += ' -use_packer_instead_of_rotamer_trials'
                        args2 += " -global_optimize"

                        prev_job_tags = [ 'REGION_0_0' ]
                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

                    else:
                        sub_job_tag = "START_FROM_START_PDB"

                        args2 = args
                        args2 += ' -s1 ' + start_pdbs[n]
                        args2 += ' -input_res1 '
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )

                        args2 += ' -use_packer_instead_of_rotamer_trials'

                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

        else:

            if TEMPLATE:

                for n in range( len( template_pdbs ) ):
                    # Good to just build the whole thing off template.
                    template_pdb = template_pdbs[ n ]
                    template_mapping = template_mappings[ template_pdb ]

                    if ( not template_continuous( i,j,template_mapping ) ): continue

                    sub_job_tag = "START_FROM_TEMPLATE_%d" % n

                    args2 = args
                    args2 += ' -s1 ' + template_pdb
                    args2 += ' -input_res1 '
                    args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )
                    args2 += ' -slice_res1 '
                    args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )
                    args2 += ' -backbone_only1'
                    args2 += ' -use_packer_instead_of_rotamer_trials'

                    setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                         fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

            if FRAGMENT_LIBRARY:
                for n in range( len( frag_files ) ):
                    # We might be at a length where we match a fragment...
                    frag_length = frag_lengths[ n ]
                    frag_file = frag_files[ n ]

                    if not ( frag_length  == L ): continue

                    sub_job_tag = "START_FROM_FRAGMENT_LIBRARY_%dMER" % frag_length

                    args2 = args
                    args2 += ' -in:file:frag_files ' + frag_file
                    args2 += ' -use_packer_instead_of_rotamer_trials'

                    args2 += ' -sample_res'
                    args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )

                    setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                         fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

            if SWA_FRAGS and (outfile_cluster in silent_files_in ) and ( L in swa_frag_lengths ) :

                sub_job_tag = "SWA_REPACK_%dMER" % L

                args2 = "%s -silent1 %s/%s" % ( args, swa_silent_file_dir, outfile_cluster )
                args2 += ' -use_packer_instead_of_rotamer_trials'
                args2 += ' -input_res1 '
                args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )
                setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                     fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



            if (L <= 2) and DENOVO: # This happens for two-residue fragments.

                sub_job_tag = "START_FROM_SCRATCH"

                args2 = args
                args2 += ' -sample_res %d %d ' % (i,j)

                setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                     fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)


        ########################################################
        # Chain closure #1 -- assume j, j+2 are sampled already,
        #  then solve for j,j+1,j+2. Also (important!), this
        #  assumes that the backbone outside the loop is fixed.
        # This is meant for O(N) loop modeling -- build forward,
        #  build backward and meet in the middle.
        ########################################################

        #if  START_FROM_PDB and loop_close and (not no_fixed_res) and ( j >= loop_start-1 and i <= loop_end+1) and not loop_force_Nsquared:
        if  START_FROM_PDB and loop_close and (not no_fixed_res) and ( j >= loop_start-1 and i <= loop_end+1):

            ###############################################################################################
            # Kinematic loop closure -- close gaps with 3 bridge residues -- might deprecate this soon
            ###############################################################################################
            if DO_KIC:
                i_prev = i+2
                j_prev = j-1
                # One parent silent file exists with exactly 3 residues to be closed. Will occur in N^2 modeling or  or if there are only 3 residues in loop [as a boundary case in O(N)],
                job_tag = "REGION_%d_%d" % ( i_prev, j_prev )
                if job_tag in all_job_tags:
                    sub_job_tag = 'START_FROM_%s_CLOSE_LOOP_KIC' % ( job_tag )

                    decoy_tag = 'S_$(Process)'

                    infile1 =  job_tag.lower()+"_sample.cluster.out"
                    args2 = '%s  -silent1 %s -tags1 %s' % (args, infile1, decoy_tag )
                    args2 += " -input_res1 "
                    args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i_prev, j_prev+1, NRES, start_res_for_input_pdb ) )

                    args2 += " -bridge_res %d %d %d" % (j,j+1,j+2)
                    args2 = add_cutpoint_closed( args2, j+1 )
                    args2 += " -global_optimize"

                    prev_job_tags = [ job_tag ]
                    setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                         fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)


                # KIC -- combine loop models from both N-terminus and C-terminus: close the loop. Again, might be deprecated soon.
                i_prev = i+2
                j_prev = j-1
                job_tag1 = "REGION_%d_%d" % ( i_prev,     loop_start-1 )
                job_tag2 = "REGION_%d_%d" % ( loop_end+1, j_prev       )

                if ( job_tag1 in all_job_tags ) and ( job_tag2 in all_job_tags) and ( j_prev >= loop_start and i_prev <= loop_end ):

                    sub_job_tag = 'START_FROM_%s_%s_CLOSE_LOOP_KIC' % ( job_tag1, job_tag2 )

                    decoy_tag = 'S_$(Process)'

                    infile1 =  job_tag1.lower()+"_sample.cluster.out"
                    args2 = '%s  -silent1 %s -tags1 %s' % (args, infile1, decoy_tag )
                    args2 += " -input_res1 "
                    args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i_prev, loop_start, NRES, start_res_for_input_pdb ) )

                    infile2 =  job_tag2.lower()+"_sample.cluster.out"
                    args2 += ' -silent2 %s ' % ( infile2 )
                    args2 += " -input_res2 "
                    args2 += make_tag_with_dashes( wrap_range_and_add_satellite(loop_end+1, j_prev+1, NRES, start_res_for_input_pdb ) )

                    args2 += " -bridge_res %d %d %d" % (j,j+1,j+2)
                    args2 = add_cutpoint_closed( args2, j+1 )

                    args2 += " -global_optimize"

                    prev_job_tags = [ job_tag1, job_tag2 ]
                    setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                         fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)


            ######################################################################
            # CCD loop closure -- close gaps with 2, 1, or 0 bridge residues
            ######################################################################
            # combine loop models from both N-terminus and C-terminus: close the loop
            if DO_CCD:
                j_prev = j
                for i_prev in [ i, i+1, i+2 ]:
                    if ( j_prev < loop_start or i_prev > loop_end ): continue

                    job_tag1 = "REGION_%d_%d" % ( i_prev,     loop_start-1 )
                    job_tag2 = "REGION_%d_%d" % ( loop_end+1, j_prev       )

                    if ( job_tag1 in all_job_tags ) and ( job_tag2 in all_job_tags ):

                        sub_job_tag = 'START_FROM_%s_%s_CLOSE_LOOP_CCD' % ( job_tag1, job_tag2 )

                        decoy_tag = 'S_$(Process)'

                        infile1 =  job_tag1.lower()+"_sample.cluster.out"
                        args2 = '%s  -silent1 %s -tags1 %s' % (args, infile1, decoy_tag )
                        args2 += " -input_res1 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i_prev, loop_start, NRES, start_res_for_input_pdb ) )

                        infile2 =  job_tag2.lower()+"_sample.cluster.out"
                        args2 += ' -silent2 %s ' % ( infile2 )
                        args2 += " -input_res2 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(loop_end+1, j_prev+1, NRES, start_res_for_input_pdb ) )

                        if ( i_prev > i ):  args2 += " -bridge_res " + make_tag( range( i, i_prev ) )
                        args2 = add_cutpoint_closed( args2, j )
                        #args2 += " -ccd_close_res %d" % j
                        args2 += " -ccd_close"
                        args2 += " -global_optimize"

                        prev_job_tags = [ job_tag1, job_tag2 ]
                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)


                # build residues and close -- again, this might be combinable with stuff above.
                # in this case, build one residue from N-terminus, sample it, and CCD close over intervening residues.
                j_prev = j-1
                for i_prev in [ i, i+1, i+2 ]:
                    job_tag = "REGION_%d_%d" % ( i_prev, j_prev )
                    if disable_sampling_of_loop_takeoff and ( j==loop_start-1 or i_prev==loop_end+1 ): continue
                    if check_in_another_loop( i_prev, loop_start, all_loop_res ): continue # special case in multiloops
                    if check_in_another_loop( j_prev, loop_start, all_loop_res ): continue # special case in multiloops
                    if job_tag in all_job_tags:
                        sub_job_tag = 'START_FROM_%s_CLOSE_LOOP_CCD' % ( job_tag )

                        decoy_tag = 'S_$(Process)'

                        infile1 =  job_tag.lower()+"_sample.cluster.out"
                        args2 = '%s  -silent1 %s -tags1 %s' % (args, infile1, decoy_tag )
                        args2 += " -input_res1 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i_prev, j_prev+1, NRES, start_res_for_input_pdb ) )

                        args2 += " -sample_res %d" % ( j )
                        if ( i_prev > i ):  args2 += " -bridge_res " + make_tag( range( i, i_prev ) )
                        args2 = add_cutpoint_closed( args2, j )
                        #args2 += " -ccd_close_res %d" % j
                        args2 += " -ccd_close"

                        args2 += " -global_optimize"

                        prev_job_tags = [ job_tag ]
                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)


                # in this case, build one residue from C-terminus, sample it, and CCD close over intervening residues.
                j_prev = j
                for i_prev in [ i+1, i+2, i+3 ]:
                    job_tag = "REGION_%d_%d" % ( i_prev, j_prev )
                    cutpoint_at_Cterm = i_prev-2
                    if disable_sampling_of_loop_takeoff and ( j ==loop_start-1 or cutpoint_at_Cterm==loop_end ): continue
                    if check_in_another_loop( i_prev, loop_start, all_loop_res ): continue # special case in multiloops
                    if check_in_another_loop( j_prev, loop_start, all_loop_res ): continue # special case in multiloops
                    if job_tag in all_job_tags:
                        sub_job_tag = 'START_FROM_%s_CLOSE_LOOP_CCD' % ( job_tag )

                        decoy_tag = 'S_$(Process)'

                        infile1 =  job_tag.lower()+"_sample.cluster.out"
                        args2 = '%s  -silent1 %s -tags1 %s' % (args, infile1, decoy_tag )
                        args2 += " -input_res1 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i_prev, j_prev+1, NRES, start_res_for_input_pdb ) )

                        args2 += " -sample_res %d" % (i_prev-1)
                        if ( i_prev > i+1 ):  args2 += " -bridge_res " + make_tag( range( i, i_prev-1 ) )
                        args2 = add_cutpoint_closed( args2, cutpoint_at_Cterm )
                        #args2 += " -ccd_close_res %d" % cutpoint_at_Cterm
                        args2 += " -ccd_close"

                        args2 += " -global_optimize"

                        prev_job_tags = [ job_tag ]
                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



        ###########################################################
        # APPEND OR PREPEND TO PREVIOUS PDB
        #  [I wonder if this could just be unified with above?]
        #
        ###########################################################
        for start_region in start_regions:

            i_prev = start_region[0]
            j_prev = start_region[1]
            prev_job_tags = [ 'REGION_%d_%d' % (i_prev,j_prev) ]
            infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)

            input_res1 =  wrap_range_and_add_satellite(i_prev, j_prev+1, NRES, start_res_for_input_pdb )

            if ( loop_close and START_FROM_PDB and DO_KIC):
                # FOLLOWING HAS NOT BEEN UPDATED AFTER IMPROVEMENTS TO
                #  LOOP CLOSURE ABOVE! -- NEED TO DOUBLE-CHECK IT.
                # Is it even necessary?
                if False:
                    if new_loop_close:
                        if ( j == j_prev ) and (i_prev == (j + min_loop_gap) ): # close the loop
                            #  sample j, j+1;  bridge_res j+2, j+3, j+4.  Note that  j+4 was previously sampled, so we can hopefully
                            #  believe its psi, omega.
                            sub_job_tag = 'START_FROM_REGION_%d_%d_CLOSE_LOOP_KIC' % ( i_prev, j_prev )

                            decoy_tag = 'S_$(Process)'

                            args2 = '%s  -silent1 %s -tags1 %s' % (args, infile, decoy_tag )
                            args2 += " -input_res1 "
                            for m in input_res1: args2 += ' %d' % m

                            args2 += " -sample_res"
                            if ( j not in fixed_res ): args2 += " %d" % j
                            args2 +=" %d" % (j+1)
                            args2 += " -bridge_res %d %d %d" % (j+2,j+3,j+4)
                            args2 = add_cutpoint_closed( args2, j+1 )

                            setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                                 fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

                    else:
                        if ( j == j_prev ) and  (i_prev == j+2) and ( j >= loop_start) and (j+2 <= loop_end):
                            #  sample j, j+1;  bridge_res j+2, j+3, j+4.  Note that  j+4 was previously sampled, so we can hopefully
                            #  believe its psi, omega.
                            sub_job_tag = 'START_FROM_REGION_%d_%d_CLOSE_LOOP_KIC' % ( i_prev, j_prev )

                            decoy_tag = '' #Should be super fast.

                            args2 = '%s  -silent1 %s ' % (args, infile )
                            args2 += " -input_res1 "
                            for m in input_res1: args2 += ' %d' % m

                            #args2 += " -sample_res"
                            #for m in wrap_range( loop_end, loop_start+1): args2 += " %d" % m

                            args2 += " -global_optimize"

                            args2 += " -bridge_res %d %d %d" % (j,j+1,j+2)
                            args2 = add_cutpoint_closed( args2, j+1 )

                            setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                                 fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

            else:
                if TEMPLATE:
                    for n in range( len( template_pdbs ) ):
                        template_pdb = template_pdbs[ n ]
                        template_mapping = template_mappings[ template_pdb ]

                        input_res2 = []
                        for m in wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ):
                            if m not in input_res1:
                                input_res2.append( m )
                        if ( not template_continuous( input_res2[0],input_res2[-1],template_mapping ) ): continue


                        sub_job_tag = 'START_FROM_REGION_%d_%d_TEMPLATE_%d' % ( i_prev, j_prev, n )

                        args2 = "%s  -silent1 %s " % (args, infile )
                        args2 += " -input_res1 "
                        for m in input_res1: args2 += ' %d' % m

                        args2 += " -s2 %s" % template_pdb

                        args2 += " -input_res2 "
                        for k in input_res2: args2 += ' %d' % k

                        args2 += " -slice_res2 "
                        for k in input_res2: args2 += ' %d' % template_mapping[ k ]

                        args2 += ' -backbone_only2'

                        args2 += ' -sample_res '
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ) )


                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)

                if FRAGMENT_LIBRARY:
                    for frag_length in frag_lengths:

                        ( startpos, endpos, num_overlap_residues ) = check_frag_overlap( i, j, i_prev, j_prev, frag_length )

                        if ( num_overlap_residues < 0 or num_overlap_residues > MAX_FRAGMENT_OVERLAP ): continue

                        sub_job_tag = 'START_FROM_REGION_%d_%d_FRAGMENT_LIBRARY_%dMER' % ( i_prev, j_prev, frag_length )

                        decoy_tag = 'S_$(Process)'

                        args2 = '%s  -silent1 %s -tags1 %s' % (args, infile, decoy_tag )
                        args2 += " -input_res1 "
                        for m in input_res1: args2 += ' %d' % m

                        frag_file = frag_files[  frag_lengths.index( frag_length) ]
                        args2 += " -in:file:frag_files %s" % frag_file

                        args2 += ' -sample_res '
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(startpos, endpos+1, NRES, start_res_for_input_pdb ) )

                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



                if SWA_FRAGS:

                    for frag_length in swa_frag_lengths:

                        ( startpos, endpos, num_overlap_residues ) = check_frag_overlap( i, j, i_prev, j_prev, frag_length )
                        if ( num_overlap_residues < 0 or num_overlap_residues > MAX_FRAGMENT_OVERLAP ): continue

                        infile_swa = "region_%d_%d_sample.cluster.out" % (startpos, endpos)
                        if infile_swa not in silent_files_in: continue

                        sub_job_tag = 'START_FROM_REGION_%d_%d_SWA_FRAGMENTS_%dMER' % ( i_prev, j_prev, frag_length )

                        decoy_tag = 'S_$(Process)'

                        args2 = '%s  -silent1 %s -tags1 %s' % (args, infile, decoy_tag )
                        args2 += " -input_res1 "
                        for m in input_res1: args2 += ' %d' % m

                        args2 += " -silent2 %s/%s" % (swa_silent_file_dir, infile_swa )
                        args2 += " -input_res2 "
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(startpos, endpos+1, NRES, start_res_for_input_pdb ) )

                        args2 += ' -sample_res '
                        args2 += make_tag_with_dashes( wrap_range_and_add_satellite(startpos, endpos+1, NRES, start_res_for_input_pdb ) )

                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



                # One final possibility. Addition of pre-existing rigid chunk (a "starting pdb" )
                if START_FROM_PDB:
                    for n in range( len(start_res_for_input_pdb) ):

                        start_res = reorder_to_be_contiguous( start_res_for_input_pdb[ n ] )

                        boundary_res = 0
                        if ( (start_res[0] == i) and (start_res[-1] in res_to_be_modeled) and  (i_prev == start_res[-1] + 1) and (j_prev == j) ):
                            boundary_res = i_prev

                        if ( ( start_res[0] in res_to_be_modeled) and (start_res[-1] == j) and (i_prev == i) and (j_prev == start_res[0]-1 ) ):
                            boundary_res = j_prev

                        if ( boundary_res == 0 ): continue

                        sub_job_tag = 'START_FROM_REGION_%d_%d_WITH_START_PDB_%d' % ( i_prev, j_prev, n )

                        args2 = "%s  -silent1 %s " % (args, infile )
                        args2 += " -input_res1 "
                        for m in input_res1: args2 += ' %d' % m

                        args2 += " -s2 %s" % start_pdbs[n]

                        input_res2 = []
                        for m in wrap_range_and_add_satellite(i, j+1, NRES, start_res_for_input_pdb ):
                            if m not in input_res1:
                                input_res2.append( m )

                        args2 += " -input_res2 "
                        for k in input_res2: args2 += ' %d' % k

                        args2 += ' -sample_res '
                        args2 += ' %d' % boundary_res

                        setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, '', \
                                                             fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)





        do_denovo = DENOVO
        if ( not DENOVO ) and len( combine_files ) == 0 and (not loop_close) and (not override):
            print "No template_jobs for ", overall_job_tag,
            print " ... rescuing the run with denovo!"
            do_denovo = 1
            max_res_to_add_denovo = 1

        for start_region in start_regions:

            i_prev = start_region[0]
            j_prev = start_region[1]
            prev_job_tags = [ 'REGION_%d_%d' % (i_prev,j_prev) ]
            infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)

            if ( abs(i - i_prev ) <= max_res_to_add_denovo and \
                 abs(j - j_prev ) <= max_res_to_add_denovo and
                 ( not loop_close ) and  do_denovo ) :

                sub_job_tag = 'START_FROM_REGION_%d_%d_DENOVO' % ( i_prev, j_prev )

                decoy_tag = 'S_$(Process)'
                args2 = '%s  -silent1 %s -tags1 %s' % (args, infile, decoy_tag )

                args2 += ' -input_res1 '
                args2 += make_tag_with_dashes( wrap_range_and_add_satellite(i_prev, j_prev+1, NRES, start_res_for_input_pdb ) )

                sample_res = []
                if ( i < i_prev and j == j_prev):
                    for m in range(i,i_prev+1):
                        #args2 += ' %d' % m
                        if (not disable_sampling_of_loop_takeoff) or (m not in fixed_res): sample_res.append( m )
                elif ( i == i_prev and j > j_prev ):
                    for m in range(j_prev,j+1):
                        #args2 += ' %d' % m
                        if (not disable_sampling_of_loop_takeoff) or (m not in fixed_res): sample_res.append( m )
                else:
                    for m in [i,j]:
                        if m not in fixed_res: sample_res.append( m )

                args2 += ' -sample_res ' + make_tag(sample_res)

                # if only sampling one residue, go ahead and sample more backbone stuff.
                if len( sample_res ) == 1: args2 = args2.replace( '-n_sample 18', '-n_sample 36' )

                setup_dirs_and_condor_file_and_tags( overall_job_tag, sub_job_tag, prev_job_tags, args2, decoy_tag, \
                                                     fid_dag, job_tags, all_job_tags, jobs_done, real_compute_job_tags, combine_files)



        print len( combine_files )

        if len( combine_files ) == 0:
            if loop_close:
                print "no loop closure allowed here. OK."
                continue
            if ( START_FROM_PDB ): continue
            print "PROBLEM: template_jobs for ", overall_job_tag
            print " Possible solution: in fragment library run, specify -min_length (the smallest region modeled) to be the frag size"
            if override:
                print "OK, override ... no exit."
                continue
            exit( 0 )

        #print combine_files

        # OUTPUT DIRECTORY
        if not( exists( overall_job_tag) ):
            system( 'mkdir -p ' + overall_job_tag )

        ################################################################
        # CLUSTER! And keep a small number of representatives (400)
        ################################################################

        outfile_cluster = overall_job_tag.lower()+'_sample.cluster.out'
        args_cluster = ' -cluster_test -silent_read_through_errors -in:file:silent %s  -in:file:silent_struct_type binary  -database %s  %s -out:file:silent %s -nstruct %d %s -score_diff_cut %8.3f' % (string.join( combine_files ), DB,  cluster_tag, outfile_cluster, FINAL_NUMBER, cluster_by_all_atom_rmsd_tag, score_diff_cut )

        if FIX_CALC_RMS_TAG:
            args_cluster += ' -working_res' + make_tag_with_dashes( res_to_be_modeled )

        condor_submit_cluster_file = 'CONDOR/REGION_%d_%d/cluster.condor' % (i,j)

        make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1, "scheduler" )

        fid_dag.write('\nJOB %s %s\n' % (overall_job_tag,condor_submit_cluster_file) )
        fid_dag.write('PARENT %s CHILD %s\n' % (string.join(job_tags),overall_job_tag) )
        fid_dag.write('SCRIPT POST %s %s %s %s\n' % (overall_job_tag, POST_PROCESS_CLUSTER_SCRIPT, outfile_cluster, overall_job_tag ) )

        all_job_tags.append(  overall_job_tag )



#####################################################################################
final_outfile = "region_FINAL.out"
if not exists( final_outfile ) and ( MIN_RES == 1 and MAX_RES == NRES ):

    last_outfiles = []
    last_jobs = []
    for i in range(1, NRES+1):
        if (i == 1): j = NRES
        else: j = i - 1
        job_tag = 'REGION_%d_%d' % (i,j)
        if job_tag in all_job_tags:
            last_jobs.append( job_tag )
            last_outfiles.append( job_tag.lower()+'_sample.cluster.out' )

    assert( len(last_outfiles) > 0 )

    print 'REGION_FINAL'
    args_cluster = ' -cluster_test -in:file:silent %s  -in:file:silent_struct_type binary  -database %s  %s -out:file:silent %s  %s -score_diff_cut %8.3f -silent_read_through_errors  -nstruct %d ' % (string.join( last_outfiles ), DB,  cluster_tag, final_outfile, cluster_by_all_atom_rmsd_tag, 2 * score_diff_cut, 10000 )

    if FIX_CALC_RMS_TAG:
        args_cluster += ' -working_res '
        args_cluster += make_tag_with_dashes( range(1,NRES+1) )

    condor_submit_cluster_file = 'CONDOR/REGION_FINAL_cluster.condor'
    make_condor_submit_file( condor_submit_cluster_file, args_cluster, 1 )

    final_job_tag = "REGION_FINAL"
    fid_dag.write('\nJOB %s %s\n' % ( final_job_tag,condor_submit_cluster_file) )
    for prev_job_tag in last_jobs:
        if ( prev_job_tag not in jobs_done ):
            fid_dag.write('PARENT %s  CHILD %s\n' % (prev_job_tag, final_job_tag) )



print
print "Total number of jobs to run (not counting clustering):", len( real_compute_job_tags )
print "Total number of final outfiles (and clustering jobs):", len( all_job_tags )




