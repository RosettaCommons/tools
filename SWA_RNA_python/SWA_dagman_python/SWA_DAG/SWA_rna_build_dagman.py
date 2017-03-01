#!/usr/bin/env python

######################################################################
from SWA_rna_build_dagman_util import *
######################################################################

from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser import SWA_parse_rosetta_options
######################################################################

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
######################################################################


#####################Parse ARGV###############################################
START_argv=copy.deepcopy(argv)

print argv
print "print len(argv)=", len(argv)

# Define sequence
fasta_file = parse_options( argv, "fasta", "" )

if (fasta_file==""): error_exit_with_message("fasta_file==\"\"")

assert( exists( fasta_file ) )
sequence = open( fasta_file  ).readlines()[1][:-1]

native_pdb = parse_options( argv, "native", "" )
input_pdbs = parse_options( argv, "s", [""] ) #These are the already build residues...such as WC basepair...
input_silent_file_list = parse_options( argv, "input_silent_file", [""] ) #Substitute a input_pdb with a silent_file..for now assume that it is the first input_pdb
input_res_full = parse_seq_num_list_option( argv, "input_res", sort_list=False )
clusterer_num_pose_kept = parse_options( argv, "clusterer_num_pose_kept", 1000 )
cutpoints_open = parse_options( argv, "cutpoint_open", [ -1 ] )
input_pathway = parse_options( argv, "input_pathway", "" )
pathway_file = parse_options( argv, "pathway_file", "" )
pathway_bulge_res = parse_options( argv, "pathway_bulge_res", [ -1 ] )
allow_normal_move_at_pathway_bulge_res = parse_options( argv, "allow_normal_move_at_pathway_bulge_res", "False"  )
allow_bulge_move_at_non_pathway_bulge_res = parse_options( argv, "allow_bulge_move_at_non_pathway_bulge_res", "False"  )
build_from_helix_to_helix_mode = parse_options( argv, "build_from_helix_to_helix_mode", "True"  )
# TODO (sripakpa): add detail Prevent cycles.

allow_bulge = parse_options( argv, "allow_bulge", "True" ) #Python boolean #change to True on April 9th, 2011
starting_elements = parse_options( argv, "starting_elements", [ -1 ]  )
allow_build_from_scratch = parse_options( argv, "allow_build_from_scratch", "False"  )
allow_bulge_right_next_to_input_helix = parse_options( argv, "allow_bulge_right_next_to_input_helix", "True"  ) #Change to true on Nov 12, 2011
dinucleotide_at_single_element_cc= parse_options( argv, "dinucleotide_at_single_element_cc", "all" )
floating_base=parse_options( argv, "floating_base", "True" ) #change to true on April 9th 2011

analytic_etable_evaluation = parse_options( argv, "analytic_etable_evaluation", "True" )

OLLM_allow_previous_clash = parse_options (argv, "OLLM_allow_previous_clash", "False")
optimize_long_loop_mode = parse_options (argv, "optimize_long_loop_mode", "False")
# optimize_long_loop_mode:
#   This is used in Sripakdeevong et al. 2011 PNAS
#   Built long loop in potentially O(n) instead of O(n**2) steps.
#   Consider the following 6 nucleotides loop:
#           N-A1-A2-A3-U4-U5-U6-N
#
#   (1) Build regions such as:
#           N-A1-              -N
#           N-A1-A2            -N
#           N-               U6-N
#           N-            U5-U6-N
#
#       But don't built regions involving both sides such as:
#           N-A1-            U6-N
#
#       The idea is that nucleotides A1 and U6 are far apart and are therefore
#       unlikely to interact with one other, so there is no point in building
#       them simultaneously.
#
#   (2) Activates combine_long_loop mode:
#
#       Consider the following region:
#           N-A1-A2       U5-U6-N
#
#       In this case A2 and U5 are only 3 nucleotides apart and perhaps
#       they are close enough to be to contact/interact. This activates a mode
#       call "combine_long_loop". In this mode, conformations from the
#       following two regions:
#           N-A1-A2            -N
#                   and
#           N-            U5-U6-N
#
#       are pass in the filterer (-algorithm filter_combine_long_loop). Only
#       conformation-pairs in which A2 and U5 are in contact will pass the
#       filter ("combine_long_loop_contact"). These conformations are then
#       output into a silent_file that will then proceed to be built "locally"
#       in a O(n**2) fashion:
#
#           N-A1-A2       U5-U6-N   (output from filterer)
#                    |
#                    V
#           N-A1-A2    U4-U5-U6-N   (intermediate subsequent sampler step)
#           N-A1-A2-U3    U5-U6-N   (intermediate subsequent sampler step)
#                    |
#                    V
#           N-A1-A2-U3-U4-U5-U6-N

OLLM_chain_closure_only = parse_options (argv, "OLLM_chain_closure_only", "False")
# OLLM_chain_closure_only :
#   This is used in Sripakdeevong et al. 2011 PNAS
#   OLLM is abbreviation for optimize_long_loop_mode
#   A strict extension of optimize_long_loop_mode (strictly run in O(n) steps
#   rather than O(n**2) steps where n is the number of loop nucleotides).
#   Should be run only in
#   conjunction with optimize_long_loop_mode.
#
#   (1) Modifies the behavior of the combine_long_loop mode:
#
#       Only regions combine regions that are 1 nt apart are combined, e.g.:
#           N-A1-A2-A3         -N
#                   and
#           N-            U5-U6-N
#
#       The structure outputted from the filterer will be:
#           N-A1-A2-A3    U5-U6-N
#
#       The chain can then be intermediately be closed in the subseqent
#       sampler step.


pure_append_prepend_only=  parse_options (argv, "pure_append_prepend_only", "False")
enforce_path_base_pairs= parse_options (argv, "enforce_path_base_pairs", [""]) #these are elements
kill_paths= parse_options(argv, "kill_paths", [""]) #these are elements

max_duplex_frame_slip= parse_options(argv, "max_duplex_frame_slip", 0)
duplex_frame_BP= parse_options(argv, "duplex_frame_BP", [ 0 ]) # these are elements

Is_valid_value_dinucleotide_at_single_element_cc(dinucleotide_at_single_element_cc)

#sample_virt_ribose_in_sep_DAG= parse_options(argv, "sample_virt_ribose_in_sep_DAG", "True") #Change this to True on Jan 29, 2012!
sample_virt_ribose_in_sep_DAG= parse_options(argv, "sample_virt_ribose_in_sep_DAG",  floating_base )
# TODO (sripakpa): add details virt_ribose sampling should not be required if no floating base moves.
# Allow better parallelization of sampler step (prevent bottleneck)
# "TODO (sripakpa): explain why have to sample virt ribose for each seperate build step that uses a particular outfile.


final_rebuild_bulge_step_only = parse_options(argv, "final_rebuild_bulge_step_only", "False") #New option Oct 22, 2011.
BMRB_chemical_shift_file = parse_options(argv, "BMRB_chemical_shift_file", "") #New option Oct 23, 2011.

allow_combine_DS_regions = parse_options(argv, "allow_combine_DS_regions", "False")
# allow_combine_DS_regions :
#   DS is abbreviation for 'double-stranded'
#   This option provides additional building steps when building double-
#   stranded-motifs (e.g. activates combine_DS_regions mode.)
#
#   Consider the following 8 nucleotides double-strand region.
#
#   N-N-A1-A2-A3-A4-N-N
#   N-N-G8-G7-G6-G5-N-N
#
#   Regular SWA building of double-strand motifs will generate partially
#   build regions such as:
#
#   N-N-A1-A2
#   N-N-G8-G7
#           and
#             A3-A4-N-N
#             G6-G5-N-N
#
#  In combine_DS_regions mode provides building sampling steps to combine the
#  two regions above into the full-length structure.


if (len(duplex_frame_BP)!=0 and len(duplex_frame_BP)!=2):
	error_exit_with_message("len(duplex_frame_BP)!=0 and len(duplex_frame_BP)!=2, duplex_frame_BP=%s " %(list_to_string(duplex_frame_BP)) )

if (len(pathway_bulge_res)!=0 and allow_bulge==False ): error_exit_with_message("len(pathway_bulge_res)!=0 but allow_bulge==False")


enforce_path_base_pair_list=[]
for base_pair_string in enforce_path_base_pairs:
	if (base_pair_string==""): continue
	base_pair=base_pair_string.split('-')
	for n in range(len(base_pair)):
		base_pair[n]=int(base_pair[n])
	if (len(base_pair)!=2): error_exit_with_message("len(base_pair (%s) )!=2" %(base_pair_string) )
	if ((base_pair[0]<base_pair[1])==False): error_exit_with_message("(base_pair[0]<base_pair[1])==False, base_pair_string=%s" %(base_pair_string) )

	#print "base_pair= ", base_pair
	enforce_path_base_pair_list.append(base_pair)
print "enforce_path_base_pair_list= ", enforce_path_base_pair_list


##################################################################
post_process_filtered_nstruct=25*clusterer_num_pose_kept #Keep enough for clustering

chemical_shift_run=True #Jan 23, 2012 Temporary/HACKY!

if (chemical_shift_run==False):

	full_length_clusterer_num_pose_kept=clusterer_num_pose_kept

	full_length_post_process_filtered_nstruct=25*full_length_clusterer_num_pose_kept  #Keep enough for clustering

	final_clusterer_num_pose_kept=10*clusterer_num_pose_kept #Want to keep the all the structures from the various full-length regions.

else: #Chemical shift run require keeping higher energies full-length structures

	#Since full length regions, not much computational cost to keep more structures.
	#Feb 13, 2012: Cap at 10000 due to memory issue.
	full_length_clusterer_num_pose_kept=min(10000, 10*clusterer_num_pose_kept)

	full_length_post_process_filtered_nstruct=int(7.5*full_length_clusterer_num_pose_kept)
	#Could be 25*full_length_clusterer_num_pose_kept, but the problem is not there is not much structure from the sampler, which only keep a maximum of 108 from each job
	#So a maximum of 108*clusterer_num_pose_kept=10.8*full_length_clusterer_num_pose_kept. I don't want to include the 25% worst percentile since these probably include poor energy structures which are probably messed up.

	final_clusterer_num_pose_kept=full_length_clusterer_num_pose_kept #Want to keep more but 10,000 structures is probably already a sizable fraction of the total memory!

print "------------------------------------------------------------------------------------"
print "clusterer_num_pose_kept=%d | full_length_clusterer_num_pose_kept=%d | final_clusterer_num_pose_kept=%d" %(clusterer_num_pose_kept, full_length_clusterer_num_pose_kept, final_clusterer_num_pose_kept)
print "post_process_filtered_nstruct=%d | full_length_post_process_filtered_nstruct=%d " %(post_process_filtered_nstruct, full_length_post_process_filtered_nstruct)
print "------------------------------------------------------------------------------------"
##################################################################


###Next Parse options that are used exclusive by the Rosetta code ####
common_args= SWA_parse_rosetta_options.get_rosetta_common_args_option(argv)
sampling_args= SWA_parse_rosetta_options.get_rosetta_samplerer_args_option(argv)
cluster_args= SWA_parse_rosetta_options.get_rosetta_clusterer_args_option(argv)

###After parsing there should only be one argv left which is the python script name (SWA_rna_build_dagman.py)####

if (len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

######################################################################################
if (exists("CONDOR/")): #This ensures that everything inside this folder is from the latest call to SWA_rna_build_dagman.py
	print "CONDOR/ folder already exist ...removing!"
	submit_subprocess("rm -r CONDOR/")

system( "mkdir CONDOR/" )

if (exists("COMMON_ARGS/")): #This ensures that everything inside this folder is from the latest call to SWA_rna_build_dagman.py
	print "COMMON_ARGS/ folder already exist ...removing!"
	submit_subprocess("rm -r COMMON_ARGS/")

system( "mkdir COMMON_ARGS/" )

if (exists(final_common_args_file())): error_exit_with_message("final_common_args_file (%s) already exist!" %(final_common_args_file()))

##################argvs shared by both Rosetta and the Python script##################
all_res=range(1,len(sequence)+1)
global_sample_res_list=[]
for seq_num in all_res:
	if (seq_num in input_res_full): continue
	global_sample_res_list.append(seq_num)

global_sample_res_list.sort()

global_sample_res_plus_bound=[]

for seq_num in all_res:
	add_res=False
	if (seq_num in global_sample_res_list): add_res=True
	if ((seq_num-1) in global_sample_res_list): add_res=True
	if ((seq_num+1) in global_sample_res_list): add_res=True

	if (add_res): global_sample_res_plus_bound.append(seq_num)

global_sample_res_plus_bound.sort()

print "all_res=                ", all_res
print "input_res_full=         ", input_res_full
print "global_sample_res_list= ", global_sample_res_list
print "global_sample_res_plus_bound= ", global_sample_res_plus_bound

if ( not analytic_etable_evaluation ): common_args += ' -analytic_etable_evaluation false'

common_args += ' -fasta %s  -global_sample_res_list %s ' %(fasta_file, list_to_string(global_sample_res_list) )

if (len( cutpoints_open ) > 0): common_args += ' -cutpoint_open %s ' %(list_to_string(cutpoints_open) )

sampling_args = ' -algorithm rna_sample  -database %s  %s ' %(DB, sampling_args)

if ( native_pdb!="" ):
	if (exists( native_pdb )==False): error_exit_with_message("exists( native_pdb (%s) )==False" %(native_pdb))
	sampling_args += '-native %s ' %(native_pdb)

cluster_args = ' -algorithm rna_cluster  -database %s  %s ' %(DB, cluster_args)
#TODO (sripakpa): Describe clusterer use cases (RMSD base comparision)

filterer_args = ' -algorithm filter_combine_long_loop -database %s -clusterer_num_pose_kept %d ' %(DB, clusterer_num_pose_kept)
#TODO (sripakpa): Describe filterer use cases (takes input two silent_files and output
#       combined structures that pass certain filters)

#########################################################
#pdb info/pose_info is a map#
#member variables
#pdb_info["tag"]   lower_elements.pdb if pdb, S_0 if in silent_file
#pdb_info["element"]=0
#pdb_info["residue_list"]=[1,2,3,4,5]

user_input_pdb_info_list=get_user_input_pdb_info_list(input_pdbs, input_silent_file_list, input_res_full, sequence)

#########################################################################################################

input_pdbs=[] #Will not uses this variable after this point, so EMPTY THEM
input_res_full=[]  #Will not uses this variable after this point, so EMPTY THEM


##DO TO: replace the coloring code....it is hard to understand

coloring = []  #maps seq_num to user_input_pdb ID
for i in range( len( sequence ) ): # sequence contain every residues, the pre-builts and the one that will be build. 	 # Be careful range(10) == [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
	color = -1
	for j in range( len( user_input_pdb_info_list ) ):
		input_pdb_info=user_input_pdb_info_list[j]
		if (i+1) in input_pdb_info["residue_list"]: color = j
	coloring.append( color )


# coloring   is        0, 0, -1, 1, 1, -1, 0, 0  #order by seq_num
# assigned_element is  0, 0,  1, 2, 2,  3, 0, 0
#
#  i.e., there are actually four moving elements, numbered 0, 1, 2, and 3
#   for assigned_element, its convenient to start numbering at zero because of
#   varius modulo operations coming up...
#

colors_so_far = []
num_elements = -1  #Start with -1 since want first member of assign_element to be zero.
color_to_element = {}  #rename to user_input_pdb_ID_to_element
element_to_color = {}  #rename to element_to_user_input_pdb_ID
element_definition = {} #map from assign_element to corresponding seq_list
assigned_element = []  #map from Seq_num to assigned_element
for i in range( len(sequence ) ): #goes from 0 len(sequence) -1

    if ( coloring[ i ] == -1 ):
        num_elements += 1
        assign_element = num_elements
    elif  (coloring[ i ]  not in colors_so_far ):
        num_elements += 1
        colors_so_far.append( coloring[ i ] )
        color_to_element[ coloring[ i ] ] = num_elements
        element_to_color[ num_elements ] = coloring[ i ]
        assign_element = num_elements
    else:
        assign_element = color_to_element[ coloring[ i ] ]

    if assign_element not in element_definition:
        element_definition[ assign_element ] = []
    element_definition[ assign_element ].append( i+1  )

    assigned_element.append( assign_element )

num_elements += 1  #This is compensate for the fact that num_elements starts out at -1....with this num_elements REALLY equal the number of elements in the system

element_definition=reorder_element_definition(element_definition, cutpoints_open, num_elements, len(sequence ))

print "#################################"
print "       DEFINITION OF ELEMENTS    "
print "#################################"
print element_definition
print " or: "
print assigned_element
print

for pdb_num in range(len(user_input_pdb_info_list) ):
	user_input_pdb_info_list[pdb_num]['element']=color_to_element[pdb_num]
	user_input_pdb_info_list[pdb_num]['i_prev']=color_to_element[pdb_num]
	user_input_pdb_info_list[pdb_num]['j_prev']=color_to_element[pdb_num]
output_pdb_info_list(user_input_pdb_info_list)

#print color_to_element
#print element_to_color



# Keep a list of all regions that we are building.
# Initialize with starting pdb -- already built!
all_job_tags = [] #input_pdb_tags or clusterer_job_tags
jobs_done = []  	 #input_pdb_tags or (clusterer_job_tags that were generated in a prior BIOX run).

last_jobs = []     #The clusterer_job_tags associated with full length pose
last_outfiles = [] #The clusterer_outfiles associate with full length pose


# Tags, etc. are now indexed by "element", 0, 1, 2, 3, etc...

build_from_scratch_tags = []

if (allow_build_from_scratch):
	for k in range( num_elements ) :

		if (k in element_to_color.keys()): continue

		if ( ((k-1) % num_elements)  in element_to_color.keys()): continue

		if ( ((k+1) % num_elements)  in element_to_color.keys()): continue

		input_file_tag = get_job_tag(k, k)
		all_job_tags.append(input_file_tag )
		jobs_done.append( input_file_tag )
		build_from_scratch_tags.append( input_file_tag )

for pdb_info in user_input_pdb_info_list:

	element=pdb_info["element"]
	input_pdb_tag = get_job_tag(element, element)

	#if specify starting_elements, that use the only starting_elements as the starting point
	if ( (len(starting_elements) > 0) and (element not in starting_elements) ): continue

	all_job_tags.append( input_pdb_tag )
	jobs_done.append( input_pdb_tag )

print "starting_regions: ", get_region_tags(user_input_pdb_info_list), " build_from_scratch_tags: ", build_from_scratch_tags


########################################################################################################
combine_DS_regions=False

if allow_combine_DS_regions:

	#The two modes are not compatible!
	if optimize_long_loop_mode:
		error_exit_with_message("BOTH optimize_long_loop_mode and allow_combine_DS_regions EQUALS True!")

	#Might not be compatible especially if len(starting_elements)>=1
	if len(starting_elements) > 0:
		error_exit_with_message("len(starting_elements) > 0 and allow_combine_DS_regions EQUALS True!")

	#Not yet tested if the two modes are compatible!
	if allow_build_from_scratch:
		error_exit_with_message("BOTH allow_build_from_scratch and allow_combine_DS_regions EQUALS True!")

	combine_DS_regions = len(cutpoints_open) != 0

	if combine_DS_regions:
		if len(user_input_pdb_info_list) != 2:
			error_exit_with_message("len(user_input_pdb_info_list)=(%s)!=2" %(len(user_input_pdb_info_list) ) )

print "allow_combine_DS_regions=%s" % allow_combine_DS_regions,
print "| combine_DS_regions=%s " % combine_DS_regions

###############################################################################

specified_pathway=False

num_elements_actual=num_elements

if (pathway_file!=""):
	if (len(pathway_bulge_res)!=0):error_exit_with_message("input_path!="" and pathway_file!="" ")
	specified_pathway=True
	(step_list, num_elements_actual)=import_step_list_from_file(pathway_file,num_elements)

	if (num_elements_actual!=num_elements and num_elements_actual!=num_elements-1):
		print "WARNING!: num_elements_check(%d)!=num_elements(%d) or num_elements-1" %(num_elements_actual,num_elements)
	else:
		print "Confirmed that num_elements_actual==num_elements or num_elements_actual==num_elements-1"

if combine_DS_regions:
	# Reverts back num_elements_actual, since combine_DS_regions guarantee that
	# every element is build.
	num_elements_actual=num_elements

###############################################################
# MAIN LOOP
###############################################################
fid_dag = open( "rna_build.dag", 'w' ) #Where the dag commands are ouputted.
fid_dag.write("DOT dag.dot\n")

total_samplerer_jobs = 0
num_bulge_moves=0

CREATED_COMMON_ARGS_FINAL=False #Feb 07, 2012: For consistency check!

# Order calculation based on number of elements modeled -- smaller fragments first.
for L in range( 2, num_elements+1 ): #from 2 to num_elements

	for k in range( num_elements ) :  #from 0 to num_elements-1

		i = k #in the range [0:num_elements-1]
		j = ( k + L - 1 ) % num_elements  #in the range [0:num_elements-1]

		# the case where i>j is only possible if there is a "outer helix/pdb" i.e. the first and last nucleotide in the sequence are part of the same element.
		if (i>j):
			first_seq_num=0
			last_seq_num=len(assigned_element)-1
			if (assigned_element[first_seq_num]!=assigned_element[last_seq_num]): continue

		##############
		if (satisfy_duplex_frame_slip(i,j, num_elements, max_duplex_frame_slip, duplex_frame_BP)==False): continue

		if (kill_paths!=[""]):
			kill_this_path=False
			for jj in range(len(kill_paths)):

				element_pairs=kill_paths[jj].split('-')

				if (len(element_pairs)!=2): error_exit_with_message("len(element_pairs)!=2" )

				if ( i==int(element_pairs[0]) and j==int(element_pairs[1]) ):
					kill_this_path=True

			if (kill_this_path==True):
				print "kill path at i=%d, j=%d " %(i,j)
				continue
		##############


		if (pure_append_prepend_only):
			if (i!=0 and j!=0): continue

		normal_move=True
		bulge_move=allow_bulge

		if (specified_pathway==True):
			good_pathway=False
			normal_move=False
			bulge_move=False

			for step in step_list:
				if (step['i']==i and step['j']==j):
					good_pathway=True
					normal_move=step['normal_move']
					bulge_move=step['bulge_move']
					break
			if (good_pathway==False): continue

		prefix = 'region_%d_%d_' % (i,j)

		clusterer_job_tag = get_job_tag(i,j)
		clusterer_outfile = prefix+'sample.cluster.out'

		###########################################
		# OUTPUT DIRECTORY
		outdir = get_job_tag(i,j)

		###########################################

		# DO THE JOBS
		start_regions = []  # a list of seq_num pair

		check_attach_verbose=False

		if (normal_move==True):

			# Prepend single nucleotide
			i_prev = (i + 1) % num_elements
			j_prev = j

			if Is_valid_attachment(i, j, i_prev, j_prev , num_elements,
								  cutpoints_open, all_job_tags, pathway_bulge_res,
								  element_definition, user_input_pdb_info_list,
								  optimize_long_loop_mode, OLLM_chain_closure_only,
								  allow_bulge_right_next_to_input_helix,
								  allow_normal_move_at_pathway_bulge_res,
								  allow_bulge_move_at_non_pathway_bulge_res,
								  dinucleotide_at_single_element_cc,
								  build_from_helix_to_helix_mode,
								  check_attach_verbose):
				start_regions.append([i_prev, j_prev])

			# Append single nucleotide
			i_prev = i
			j_prev = (j - 1) % num_elements

			if Is_valid_attachment(i, j, i_prev, j_prev , num_elements,
								  cutpoints_open, all_job_tags, pathway_bulge_res,
								  element_definition, user_input_pdb_info_list,
								  optimize_long_loop_mode, OLLM_chain_closure_only,
								  allow_bulge_right_next_to_input_helix,
								  allow_normal_move_at_pathway_bulge_res,
								  allow_bulge_move_at_non_pathway_bulge_res,
								  dinucleotide_at_single_element_cc,
								  build_from_helix_to_helix_mode,
								  check_attach_verbose):
				start_regions.append([i_prev, j_prev])

		# Dinucleotide (e.g. Bulge + floating base) move.
		if bulge_move:

			# Via prepend.
			i_prev = (i + 2) % num_elements
			j_prev = j

			if Is_valid_attachment(i, j, i_prev, j_prev , num_elements,
								  cutpoints_open, all_job_tags, pathway_bulge_res,
								  element_definition, user_input_pdb_info_list,
								  optimize_long_loop_mode, OLLM_chain_closure_only,
								  allow_bulge_right_next_to_input_helix,
								  allow_normal_move_at_pathway_bulge_res,
								  allow_bulge_move_at_non_pathway_bulge_res,
								  dinucleotide_at_single_element_cc,
								  build_from_helix_to_helix_mode,
								  check_attach_verbose):
				start_regions.append([i_prev, j_prev])

			# Via append.
			i_prev = i
			j_prev = (j - 2) % num_elements

			if Is_valid_attachment(i, j, i_prev, j_prev , num_elements,
								  cutpoints_open, all_job_tags, pathway_bulge_res,
								  element_definition, user_input_pdb_info_list, \
								  optimize_long_loop_mode, OLLM_chain_closure_only,
								  allow_bulge_right_next_to_input_helix, \
								  allow_normal_move_at_pathway_bulge_res,
								  allow_bulge_move_at_non_pathway_bulge_res, \
								  dinucleotide_at_single_element_cc,
								  build_from_helix_to_helix_mode,
								  check_attach_verbose):
				start_regions.append([i_prev, j_prev])

		samplerer_post_process_outfile_list = []
		samplerer_tag_list = []
		job_specific_common_args_INSTANCE=""


		#################Create sampling JOBS##################################

		for start_region in start_regions:

			i_prev = start_region[0]
			j_prev = start_region[1]

			##################################################################
			if (satisfy_enforce_path_base_pair(enforce_path_base_pair_list, i, i_prev, j, j_prev, allow_bulge, num_elements)==False): continue

			###################################################################

			# Modes: (1) regular addition of a residue; (2) combine long loop.
			for mode_num in [1,2]:

				if (optimize_long_loop_mode==False and mode_num==2): continue

				combine_long_loop=False

				if (mode_num==2):
					combine_long_loop=True
					if ( get_job_tag(0,j_prev) not in all_job_tags):
						print "combine_long_loop=True, early break starting_region: %s is not in all_job_tags " %(get_job_tag(0,j_prev))
						continue

					if ( get_job_tag(i_prev,0) not in all_job_tags):
						print "combine_long_loop=True, early break starting_region: %s is not in all_job_tags " %(get_job_tag(i_prev,0))
						continue

				parent_job_tag=get_job_tag(i_prev,j_prev)

				if (combine_long_loop):
					start_region_string= '0_%d_AND_%d_0' %(j_prev,i_prev)
				else:
					start_region_string= '%d_%d' % (i_prev, j_prev)

				if (combine_long_loop==False):
					if (parent_job_tag not in all_job_tags): continue

				######################################################
				# Start with "fixed", pre-constructed pose
				######################################################
				# Job input depends on whether we are starting from pdb, silent file, etc.

				job_specific_common_args = common_args
				job_specific_sampling_args= sampling_args
				job_specific_filterer_args= filterer_args

				if (combine_long_loop):
					pose_info_list=get_modeled_pose_info_list_combine_long_loop(i, j ,i_prev, j_prev, num_elements, user_input_pdb_info_list, element_definition, optimize_long_loop_mode, OLLM_chain_closure_only)

				else:

					pose_info_list=[get_modeled_pose_info(i, j ,i_prev, j_prev, num_elements, build_from_scratch_tags, user_input_pdb_info_list, element_definition)]

					if (optimize_long_loop_mode):

						# Prevent building from BOTH side (will be taken care
						# of by the combine_long_loop step).
						if i_prev != 0 or j_prev != 0: #not the first building step.
							if i != 0 and i_prev == 0: pose_info_list=[]
							if j != 0 and j_prev == 0: pose_info_list=[]

						# STRICTER CONDITION, ONLY pure append or pure prepend..
						if (OLLM_chain_closure_only):
							if ( i != 0 and j != 0): pose_info_list=[]

				if (len(pose_info_list)==0 ): #this is mainly for the optimize_long_loop_mode...
					print "Skipping region %s since len(pose_info_list)=0, combine_long_loop=%s" %(get_job_tag(i,j), combine_long_loop)
					continue

				######################################################
				# Then define the moving element.
				######################################################
				if ( i == i_prev ):
					moving_element = j
				elif ( j == j_prev):
					moving_element = i
				else:
					error_exit_with_message("(i != i_prev) and  (j != j_prev)")

				################################################################################################
				sample_res_list=[]
				sample_virtual_ribose_list=[] #list of possible virtual_ribose positions that would need to be rebuilt before the current jobs.

				sample_res = -1
				if ( i == i_prev ):
					sample_res_list.append(element_definition[ moving_element ][0]) # append: beginning residue of moving element
				else:
					sample_res_list.append(element_definition[ moving_element ][-1]) # prepend: last residue of moving element


				if (len( element_definition[ moving_element ] ) != 1): # A user input chunk!
					pose_info_list.append( user_input_pdb_info_list[ element_to_color[ moving_element ] ] )

					# P = prepend, A = append. Perhaps should make this more explicit in the command line.
					if ( i == i_prev ): #This is an append. We want to mark the VIRT ribose in the moving element which may have been prepended in an earlier job.
						sample_virtual_ribose_list.append("%d-%s" %(element_definition[ moving_element ][ 0], 'P' ) ) #Corrected from -1 to 0 on Sept 17, 2011
					else:
						sample_virtual_ribose_list.append("%d-%s" %(element_definition[ moving_element ][-1], 'A' ) ) #Corrected from 0 to -1 on Sept 17, 2011



				##Bulge res
				Bulge_res=False
				if ( i_prev == (i + 2) % num_elements and j_prev == j): #prepend + 1 bulge case

					sample_res_list.append(element_definition[ (i + 1) % num_elements ][0])
					Bulge_res=True

				elif ( i_prev == i and j_prev == (j - 2) % num_elements): #append + 1 bulge case

					sample_res_list.append(element_definition[ (j - 1) % num_elements ][0])
					Bulge_res=True

				else: #OK as a consistency check, check that the building step is a normal append or prepend (no bulge res)

					Bulge_res=False
					if ( (( i_prev == i and j_prev == (j - 1) % num_elements ) or ( i_prev == (i + 1) % num_elements and j_prev == j ) ) == False):
						error_exit_with_message("Invalid i(%d) , j(%d)" %(i,j) )

				job_specific_common_args += ' -sample_res %s ' %(list_to_string(sample_res_list) )

				if (Bulge_res and floating_base): job_specific_common_args += ' -floating_base true '

				#July 21, 2011 For special case of -bulge_res and -floating_base true [dinucleotide move],
				# the code has been optimized such that no speed gained parallelization from separating the VIRT ribose
				if ( not Bulge_res  or  not floating_base ): #For the other sampling cases, the parallelization does speed up the code.
					if ( i == i_prev ): # the moving element is appended. there may be a virtual residue at the 3' end of the template
						sample_virtual_ribose_list.append("%d-%s" %( element_definition[ (j_prev) % num_elements ][-1], 'A' ) )
					else:
						sample_virtual_ribose_list.append("%d-%s" %( element_definition[ (i_prev) % num_elements ][ 0], 'P' ) )

				if not is_release_mode():
					# Release_mode specifies that features have been tested and
					# is release for usage by the wider community.
					# Features such as consistency check should only be while
					# testing the code and not in the released version of the
					# code.
					if (sample_virt_ribose_in_sep_DAG): job_specific_sampling_args += ' -sampler_assert_no_virt_ribose_sampling true '

				# Special --> close cutpoint.
				if ( L == num_elements ):  #In the case, the element wrap into 1 exact cycle
					boundary_res = get_boundary_res( j, assigned_element )
					if ( boundary_res not in cutpoints_open):
						if (boundary_res==len(sequence)):
							print "boundary_res(%s)==len(sequence)(%s), cutpoint_closed not add" %(boundary_res, len(sequence) )
						else:
							job_specific_common_args += ' -cutpoint_closed %d ' % boundary_res
							sample_virtual_ribose_list.append("%d-%s" %(boundary_res,   'A' ) ) #VIRT ribose at the five_prime CB
							sample_virtual_ribose_list.append("%d-%s" %(boundary_res+1, 'P' ) ) #VIRT ribose at the three_prime CB
				################################################################################################

				filterer_mode=""

				if (combine_long_loop):

					Is_combine_long_loop_chain_closure_step=False
					if ( ( i% num_elements) == ( (j+1) % num_elements ) ): Is_combine_long_loop_chain_closure_step=True

					print "region i=%d j=%d, Is_combine_long_loop_chain_closure_step= %s" %(i,j,Is_combine_long_loop_chain_closure_step)

					####CHANGE to filter for clash instead of contact at Closure closure step even if OLLM_chain_closure_only==False Nov 13, 2010
					if OLLM_chain_closure_only or Is_combine_long_loop_chain_closure_step:
						if OLLM_allow_previous_clash:
							filterer_mode="combine_long_loop_allow_previous_clash"
						else:
							filterer_mode="combine_long_loop_clash"
					else:
						filterer_mode="combine_long_loop_contact"

					if sample_virt_ribose_in_sep_DAG:
						job_specific_filterer_args += " -filterer_undercount_sugar_rotamers true "
						# Generally filterer will select and output a maximum
						# counts of low energy conformations.
						# -filterer_undercount_ribose_rotamers options tells
						# the filterer that two conformations that differ only
						# by the coordinates of theirs virtual ribose/sugar
						# should be treated as a single count. In the mode,
						# the number of conformations actually outputted
						# by the filterer with exceed the 'maximum counts'

				job_specific_common_args=add_input_res_to_common_args(job_specific_common_args, pose_info_list)

				#########################################################################################################################
				if (job_specific_common_args_INSTANCE==""):
					job_specific_common_args_INSTANCE=job_specific_common_args
					make_common_args_file(job_specific_common_args_INSTANCE, 'common_args_' + basename(clusterer_outfile) )

				if ( (L == num_elements_actual) and (CREATED_COMMON_ARGS_FINAL==False ) ):

						if ( exists(final_common_args_file()) ): error_exit_with_message("final_common_args_file (%s) already exist!" %(final_common_args_file())) #Simple consistency check

						make_common_args_file(job_specific_common_args_INSTANCE, 'common_args_region_FINAL.out' )

						CREATED_COMMON_ARGS_FINAL=True

				#######################################The actual dag submission!#####################################################
				submit_DAG=True

				if (exists( clusterer_outfile )): #If clusterer_outfile already exist it means this job is already DONE by a prior SWA_rna_build_dagman.py call!
					submit_DAG=False

					if (clusterer_job_tag not in all_job_tags): all_job_tags.append(  clusterer_job_tag )
					if (clusterer_job_tag not in jobs_done): jobs_done.append( clusterer_job_tag )

					if ( L == num_elements ):  #The case where the element wrap into 1 exact cycle
						if (clusterer_job_tag not in last_jobs): 		  last_jobs.append( clusterer_job_tag )
						if (clusterer_outfile not in last_outfiles): last_outfiles.append( clusterer_outfile )

				job_tag = 'REGION_%d_%d_START_FROM_REGION_%s' % ( i, j, start_region_string)

				output_foldername='%s/START_FROM_REGION_%s' %(outdir , start_region_string)

				samplerer_post_process_filter_outfile = get_sampler_post_process_filter_outfile(output_foldername)

				samplerer_post_process_outfile_list.append( samplerer_post_process_filter_outfile )

				if ( exists(samplerer_post_process_filter_outfile) ): submit_DAG=False

				if (submit_DAG):

					filterer_parent_job_tag_list=[]
					samplerer_parent_job_tag_list=[]

					################################################################################################################################################################
					if ( (sample_virt_ribose_in_sep_DAG) and (len(sample_virtual_ribose_list)>0 ) ):

						num_silent_files = get_num_silent_files(pose_info_list)

						for n in range(len(pose_info_list)):

							prev_clusterer_job_tag = get_job_tag(pose_info_list[n]['i_prev'], pose_info_list[n]['j_prev'])

							# rhiju: I think this changes the name of the input pdb/input file for the current job
							# to take in the output of the virtual ribose sampler job. Pretty elegant, parin!
							(pose_info_list[n]['silent_file'], sample_virt_ribose_DAG_already_done) = create_sampled_virt_ribose_silent_file(pose_info_list[n], job_tag, prev_clusterer_job_tag, all_job_tags, jobs_done, sample_virtual_ribose_list, native_pdb, fid_dag, num_silent_files)

							if (not sample_virt_ribose_DAG_already_done):
								virt_ribose_sampler_job_tag=get_virt_ribose_sampler_job_tag(pose_info_list[n], job_tag)

								samplerer_parent_job_tag_list.append( virt_ribose_sampler_job_tag )
								filterer_parent_job_tag_list.append( virt_ribose_sampler_job_tag )

					################################################################################################################################################################

					create_samplerer_dag_job_file(job_specific_common_args, job_specific_sampling_args, job_specific_filterer_args, \
									      filterer_mode, pose_info_list, output_foldername, job_tag, prefix)


					if ( L == num_elements_actual ):
						ACT_post_process_filtered_nstruct = full_length_post_process_filtered_nstruct
					else:
						ACT_post_process_filtered_nstruct = post_process_filtered_nstruct

					# PRE and POST scripts should be defined after JOB declaration in condor DAGman file!
					#submit_sampling_pre_and_post_process(fid_dag, job_tag, pose_info_list, output_foldername, ACT_post_process_filtered_nstruct)

					if (combine_long_loop):

						#Might want to add parent_job_tag to the two parente_jog_tag list in the future for score comparison.
						filterer_parent_job_tag_list.append(get_job_tag(0,j_prev)) #prev_job_tag_append
						filterer_parent_job_tag_list.append(get_job_tag(i_prev,0)) #prev_job_tag_prepend


						samplerer_parent_job_tag_list.append(get_job_tag(0,j_prev)) #prev_job_tag_append
						samplerer_parent_job_tag_list.append(get_job_tag(i_prev,0)) #prev_job_tag_prepend

					else:

						samplerer_parent_job_tag_list.append(parent_job_tag)


					update_dag_dependency( fid_dag, job_tag, pose_info_list, filterer_parent_job_tag_list, samplerer_parent_job_tag_list, all_job_tags, jobs_done)

					submit_sampling_pre_and_post_process(fid_dag, job_tag, pose_info_list, output_foldername, ACT_post_process_filtered_nstruct)

					samplerer_tag_list.append(job_tag)

					if (Bulge_res==True): num_bulge_moves+=1

				################################################################################################

		#############################Clusterer of region_i_j#####################################

		total_samplerer_jobs += len(samplerer_tag_list)

		if ( len(samplerer_tag_list)>0 and exists(clusterer_outfile) ):
			error_exit_with_message( "len(samplerer_tag_list)>0 but clusterer_outfile (%s) already exists!" %(clusterer_outfile) )

		if ( len(samplerer_post_process_outfile_list)>0 and exists(clusterer_outfile)==False ):

			modeled_elements = get_modeled_elements( i, j , num_elements)  #This is a list that runs from i to j in the modulo num_elements sense

			print '\n Modeling elements: ', modeled_elements, '\n job_specific_common_args_INSTANCE: ', job_specific_common_args_INSTANCE

			#consistency check
			if (clusterer_outfile != prefix+'sample.cluster.out'):
				error_exit_with_message("clusterer_outfile=%s=!=%s=prefix+'sample.cluster.out" %(clusterer_outfile,prefix+'sample.cluster.out') )

			condor_submit_cluster_file = get_condor_submit_file(clusterer_job_tag+ "_cluster", "CLUSTERER")


			if ( L == num_elements_actual ):
				# ACT is shorthand for actual.
				ACT_clusterer_num_pose_kept=full_length_clusterer_num_pose_kept
			else:
				ACT_clusterer_num_pose_kept=clusterer_num_pose_kept
			#############################

			create_clusterer_dag_job_file( cluster_args, job_specific_common_args_INSTANCE, samplerer_post_process_outfile_list, \
																	  clusterer_outfile, condor_submit_cluster_file, ACT_clusterer_num_pose_kept)
			#update clusterer dag dependency
			fid_dag.write('\nJOB %s %s\n' % (clusterer_job_tag ,condor_submit_cluster_file) )
			fid_dag.write('PARENT %s CHILD %s\n' % (string.join( samplerer_tag_list ), clusterer_job_tag) )

			all_job_tags.append( clusterer_job_tag )


			if ( L == num_elements_actual ):
				last_jobs.append( clusterer_job_tag )
				last_outfiles.append( clusterer_outfile )


if (CREATED_COMMON_ARGS_FINAL==False): error_exit_with_message("CREATED_COMMON_ARGS_FINAL==False!")

#################################Clusterer of region_FINAL.out##########################################################

if (len(last_jobs)!=0):

	final_outfile = standard_region_FINAL()
	final_job_tag = standard_final_job_tag()

	if (exists(final_outfile)):

		if (final_rebuild_bulge_step_only==False):

			if ( len(all_job_tags)!=len(jobs_done) ):
				error_exit_with_message("len(all_job_tags)=(%s)!=(%s)=len(jobs_done) but final_outfile (%s) already exist!" %(len(all_job_tags),len(jobs_done), final_outfile))

			for clusterer_outfile in last_outfiles:
				if (exists(clusterer_outfile)==False): error_exit_with_message("final_outfile (%s) exist but clusterer_outfile (%s) doesn't exist!" %(final_outfile, clusterer_outfile) )

		print "WARNING: final_outfile (%s) already exists!" %(final_outfile)

	else:

		print '\nFINAL_COMMON_ARGS_LINE: ', get_FINAL_COMMON_ARGS_LINE()

		condor_submit_cluster_file = get_condor_submit_file(final_job_tag+ "_cluster", "CLUSTERER")

		create_clusterer_dag_job_file( cluster_args, get_FINAL_COMMON_ARGS_LINE(), last_outfiles, final_outfile, condor_submit_cluster_file, final_clusterer_num_pose_kept )

		fid_dag.write('\nJOB %s %s\n' % ( final_job_tag, condor_submit_cluster_file) )

		for clusterer_job_tag in last_jobs:

			if ( clusterer_job_tag not in all_job_tags): error_exit_with_message("clusterer_job_tag (%s) is not in all_job_tags" %(clusterer_job_tag) )
			if ( clusterer_job_tag not in jobs_done ): fid_dag.write('PARENT %s  CHILD %s\n' % (clusterer_job_tag, final_job_tag) )


#####################################################################################################################
print
print "Total number of samplerer_jobs to run = %d" %(total_samplerer_jobs)
print "Total number of 'other' jobs to run = %d" %(len( all_job_tags )-len(jobs_done) ) #This is not accurate!
print "Num bulge move = %d"  %(num_bulge_moves)
print
print "jobs_already_done   = ", jobs_done
print
print
print "all_job_tags        = ", all_job_tags
print

#####################################################################################################################



if combine_DS_regions:

	print_title_text("COMBINE_DOUBLE_STRAND_REGIONS")

	total_DS_samplerer_jobs=0

	DS_clusterer_job_tag_list=[]
	DS_clusterer_outfile_list=[]

	#if (len(user_input_pdb_info_list)!=2): error_exit_with_message("combine_DS_regions but len(user_input_pdb_info_list)!=2)" )

	region_list=[]
	create_input_silent_for_DS_combine_job_tags=[]

	combine_DS_folder="COMBINE_DS_REGIONS/"

	submit_subprocess("mkdir -p %s" %(combine_DS_folder))

	for clusterer_job_tags in all_job_tags:

		if (sample_virt_ribose_in_sep_DAG and clusterer_job_tags[0:20]=="VIRT_RIBOSE_SAMPLER_"): continue

		if (clusterer_job_tags.split("_")[0]!="REGION"):
			error_exit_with_message('clusterer_job_tags.split("_")[0]!="REGION", job_tags.split("_")[0]=%s ' %(clusterer_job_tags.split("_")[0]) )

		i=int(clusterer_job_tags.split("_")[1])
		j=int(clusterer_job_tags.split("_")[2])

		curr_region=[i,j]

		if (region_list.count(curr_region)!=0): error_exit_with_message("curr_region (%s) already exist in region_list" %(list_to_string(curr_region)) ) #ENSURE THERE IS NO DUPLICATES

		region_list.append(curr_region)

	print
	print "region_list= ", region_list
	print

	for lower_region in region_list:
		for upper_region in region_list:

			#if (lower_region[0]<=lower_region[1] and (lower_region[0]!=0) ): continue #one of the region must include the outer_element 0...define this as lower_region.
			#if (upper_region[0]>upper_region[1] or (lower_region[0]==0) ): continue

			######################################June 15, 2011###############################################
			if (0 not in get_modeled_elements(lower_region[0], lower_region[1], num_elements)): continue #lower element should contain the 0 element.
			if (0 in get_modeled_elements(upper_region[0], upper_region[1], num_elements)): continue     #upper element should not contain the 0 element.

			if (lower_region[0]==lower_region[1]): continue #This is already accounted in the standard sampling above. (The lower_helical stub)
			if (upper_region[0]==upper_region[1]): continue #This is already accounted in the standard sampling above. (The upper_helical stub)
			###################################################################################################

			# Focusing on chain closure
			if ( ( (lower_region[0]-1) % num_elements ) !=	( (upper_region[1]) % num_elements ) ): continue

			if ( ( (lower_region[1]+1) % num_elements ) !=	( (upper_region[0]) % num_elements ) ): continue

			print "lower_region= ", lower_region, " upper_region= ", upper_region


			# Following are explicitly disallowed above... not sure why checked again? -- rhiju.
			if (upper_region[0]==upper_region[1]):
				if (upper_region[0]==0): error_exit_with_message("upper_region[0]==0" )

				#if (user_input_pdb_info_list[1]['element']!=upper_region[0]):
				#	error_exit_with_message("user_input_pdb_info_list[1]['element']=(%s)!=(%s)=upper_region[0]" %(user_input_pdb_info_list[1]['element'], upper_region[0] ) )

			if (lower_region[0]==lower_region[1]):
				if (lower_region[0]!=0): error_exit_with_message("lower_region[0]!=0, lower_region[0]=%s" %(lower_region[0]) )

				if (user_input_pdb_info_list[0]['element']!=0):
					error_exit_with_message("user_input_pdb_info_list[0]['element']=(%s)!=(%s)=0" %(user_input_pdb_info_list[0]['element'], 0 ) )


			##treat case where lower_region is lower_element_pdb or upper_region is upper_region_pdb seperately.
			#if (lower_region[0]!=lower_region[1] and upper_region[0]!=upper_region[1]):

			clusterer_job_tag = 'COMBINE_REGION_%d_%d_and_%d_%d' % (lower_region[0],lower_region[1],upper_region[0],upper_region[1])

			combine_DS_folder=get_combine_DS_REGION_FOLDER()

			submit_subprocess("mkdir -p %s " %(combine_DS_folder) )

			outdir = '%s/%s/' % (combine_DS_folder, clusterer_job_tag)

			#if (exists( outdir )==False): submit_subprocess( 'mkdir -p '+outdir ) #March 18, 2011

			clusterer_outfile = combine_DS_folder + clusterer_job_tag.lower()+'_sample.cluster.out'

			###check if job is already done###
			if exists( clusterer_outfile ):
				all_job_tags.append( clusterer_job_tag )
				jobs_done.append( clusterer_job_tag )

				DS_clusterer_job_tag_list.append( clusterer_job_tag )
				DS_clusterer_outfile_list.append( clusterer_outfile )

				continue
			###################################
			lower_pose_info=get_modeled_pose_info(0, 0 ,lower_region[0], lower_region[1], num_elements, build_from_scratch_tags, user_input_pdb_info_list, element_definition )
			upper_pose_info=get_modeled_pose_info(0, 0 ,upper_region[0], upper_region[1], num_elements, build_from_scratch_tags, user_input_pdb_info_list, element_definition )

			cluster_silent_file_before_DS_combine=True

			if (cluster_silent_file_before_DS_combine):
				##June 13, 2011...OK might not want to use region_silent file directly since too structs combination.
				##For example if clusterer_num_pose_kept=1000, then could contain 1000*1000=1 millions structs combination.
				##So instead cluster the silent_region that will be used as input for the DS_combine sampling!
				lower_pose_info['silent_file']=	\
							create_input_silent_file_for_DS_combine(lower_region, lower_pose_info, cluster_args, all_job_tags, jobs_done, \
																										create_input_silent_for_DS_combine_job_tags, fid_dag, clusterer_num_pose_kept, global_sample_res_plus_bound)
				upper_pose_info['silent_file']=	\
							create_input_silent_file_for_DS_combine(upper_region, upper_pose_info, cluster_args, all_job_tags, jobs_done, \
																										create_input_silent_for_DS_combine_job_tags, fid_dag, clusterer_num_pose_kept, global_sample_res_plus_bound)

			###################################
			# Can close circle in two ways:
			# append  = continuous at lower-3'/upper-5' boundary, chainbreak at upper-3'/lower-5' boundary
			# prepend = continuous at upper-3'/lower-5' boundary, chainbreak at lower-3'/upper-5' boundary

			attach_side_list=["append", "prepend"]

			samplerer_post_process_outfile_list = []
			samplerer_tag_list = []

			job_specific_common_args_INSTANCE=""

			PRISTINE_lower_pose_info=copy.deepcopy(lower_pose_info)
			PRISTINE_upper_pose_info=copy.deepcopy(upper_pose_info)


			for side in attach_side_list:

				lower_pose_info=copy.deepcopy(PRISTINE_lower_pose_info)
				upper_pose_info=copy.deepcopy(PRISTINE_upper_pose_info)

				job_tag = '%s_BY_%s' % ( clusterer_job_tag, side.upper()) #samplerer job_tag

				output_foldername='%s/BY_%s' %(outdir ,  side.upper())

				job_specific_common_args = common_args
				job_specific_sampling_args= sampling_args + " -minimizer_move_map_option DS_COMBINE "
				job_specific_filterer_args= filterer_args

				if (cluster_silent_file_before_DS_combine):
					lower_input_region_job_tag=get_create_input_ds_combine_region_job_tag(lower_region)
					upper_input_region_job_tag=get_create_input_ds_combine_region_job_tag(upper_region)

				else:
					lower_input_region_job_tag=get_job_tag(lower_region[0] ,lower_region[1] )
					upper_input_region_job_tag=get_job_tag(upper_region[0] ,upper_region[1] )

				parent_job_tag_list=[lower_input_region_job_tag, upper_input_region_job_tag]

				# Sample_virt_ribose_in_sep_DAG for combine double_strand regions
				if sample_virt_ribose_in_sep_DAG:

					if not is_release_mode():
						# Features such as consistency check should only be while
						# testing the code and not in the released version of the
						# code.
						job_specific_sampling_args += ' -sampler_assert_no_virt_ribose_sampling true '

					num_silent_files=get_num_silent_files([lower_pose_info, upper_pose_info])

					if (num_silent_files!=2): error_exit_with_message("num_silent_files=(%s)!=2" %(num_silent_files))

					######################################
					sample_virtual_ribose_list=[]
					sample_virtual_ribose_list.append("%d-%s" %(element_definition[ lower_region[0] ][ 0], 'P') )
					sample_virtual_ribose_list.append("%d-%s" %(element_definition[ lower_region[1] ][-1], 'A') )

					(lower_pose_info['silent_file'], sample_virt_ribose_DAG_already_done) = create_sampled_virt_ribose_silent_file(lower_pose_info, job_tag, lower_input_region_job_tag ,   \
																			     all_job_tags, jobs_done, sample_virtual_ribose_list, native_pdb, fid_dag, num_silent_files)

					if (not sample_virt_ribose_DAG_already_done): parent_job_tag_list.append( get_virt_ribose_sampler_job_tag(lower_pose_info, job_tag) )
					######################################

					######################################
					sample_virtual_ribose_list=[]
					sample_virtual_ribose_list.append("%d-%s" %(element_definition[ upper_region[0] ][ 0], 'P') )
					sample_virtual_ribose_list.append("%d-%s" %(element_definition[ upper_region[1] ][-1], 'A') )

					(upper_pose_info['silent_file'], sample_virt_ribose_DAG_already_done)=create_sampled_virt_ribose_silent_file(upper_pose_info, job_tag, upper_input_region_job_tag ,   \
																																					  all_job_tags, jobs_done, sample_virtual_ribose_list, native_pdb, fid_dag, num_silent_files)

					if (not sample_virt_ribose_DAG_already_done): parent_job_tag_list.append( get_virt_ribose_sampler_job_tag(upper_pose_info, job_tag) )
					######################################

				###########################################################

				pose_info_list=[copy.deepcopy(lower_pose_info), copy.deepcopy(upper_pose_info)]

				#assigned_element = []  #map from Seq_num to assigned_element

				if (side=="append"):

					sample_res     = element_definition[ upper_region[0] ][0]
					cutpoint_closed = element_definition[ upper_region[1] ][-1]

				else: #prepend

					sample_res     = element_definition[ lower_region[0] ][0]
					cutpoint_closed = element_definition[ lower_region[1] ][-1]

				job_specific_common_args += ' -sample_res %s ' %(sample_res)
				job_specific_common_args += ' -cutpoint_closed %s ' %(cutpoint_closed)


				################################################################################################################################

				samplerer_post_process_filter_outfile = get_sampler_post_process_filter_outfile(output_foldername)

				samplerer_post_process_outfile_list.append( samplerer_post_process_filter_outfile )

				if ( exists(samplerer_post_process_filter_outfile)==False ):

					job_specific_common_args=add_input_res_to_common_args(job_specific_common_args, pose_info_list)

					create_samplerer_dag_job_file(job_specific_common_args, job_specific_sampling_args, job_specific_filterer_args, \
									      "combine_DS_regions", pose_info_list, output_foldername, job_tag, clusterer_job_tag.lower() + "_" )

					if ( L == num_elements_actual ):
						ACT_post_process_filtered_nstruct=full_length_post_process_filtered_nstruct
					else:
						ACT_post_process_filtered_nstruct=post_process_filtered_nstruct

					# PRE and POST scripts should be defined after JOB declaration in condor DAGman file!
					#submit_sampling_pre_and_post_process(fid_dag, job_tag, pose_info_list, output_foldername, ACT_post_process_filtered_nstruct)

					update_dag_dependency( fid_dag, job_tag, pose_info_list, parent_job_tag_list, parent_job_tag_list, all_job_tags, jobs_done)
					samplerer_tag_list.append( job_tag )

					submit_sampling_pre_and_post_process(fid_dag, job_tag, pose_info_list, output_foldername, ACT_post_process_filtered_nstruct)

					if (len(job_specific_common_args_INSTANCE)==0): job_specific_common_args_INSTANCE=job_specific_common_args
				################################################################################################################################


			#############################Clusterer of COMBINE_REGION_lower_i_lower_j_and_upper_i_upper_j#####################################
			total_DS_samplerer_jobs+=len(samplerer_tag_list)


			if ( len(samplerer_tag_list) != 0 ): #This particular region_i_j was sampled.

				print '\n Modeling REGION: ', clusterer_job_tag, '\n DS_job_specific_common_args_INSTANCE: ', job_specific_common_args_INSTANCE

				condor_submit_cluster_file = get_condor_submit_file(clusterer_job_tag+ "_cluster", "CLUSTERER")

				create_clusterer_dag_job_file( cluster_args, job_specific_common_args_INSTANCE, samplerer_post_process_outfile_list, \
																		 clusterer_outfile, condor_submit_cluster_file, full_length_clusterer_num_pose_kept)

				#update clustere dag dependency
				fid_dag.write('\nJOB %s %s\n' % (clusterer_job_tag ,condor_submit_cluster_file) )
				fid_dag.write('PARENT %s CHILD %s\n' % (string.join( samplerer_tag_list ), clusterer_job_tag) )

				all_job_tags.append( clusterer_job_tag )


				DS_clusterer_job_tag_list.append( clusterer_job_tag )
				DS_clusterer_outfile_list.append( clusterer_outfile )

	#if (total_DS_samplerer_jobs==0): error_exit_with_message("combine_DS_regions=True but total_DS_samplerer_jobs=0!")

	#################################Clusterer of combine_ds_region_FINAL.out##########################################################

	combine_DS_final_outfile = "%s/combine_ds_region_FINAL.out" %(combine_DS_folder)

	if (exists(combine_DS_final_outfile)):

		if (final_rebuild_bulge_step_only==False):

			for DS_clusterer_outfile in DS_clusterer_outfile_list:
				if (exists(DS_clusterer_outfile)==False):
					error_exit_with_message("combine_DS_final_outfile (%s) already exist BUT DS_clusterer_outfile (%s) doesn't exist!" %(combine_DS_final_outfile, DS_clusterer_outfile) )


		print "WARNING: combine_DS_final_outfile (%s) already exists!" %(combine_DS_final_outfile)

	else:

		combine_DS_final_job_tag = "COMBINE_DS_REGION_FINAL"

		print "\nCREATE %s JOB" %(combine_DS_final_job_tag)

		combine_DS_final_condor_file = get_condor_submit_file(combine_DS_final_job_tag+ "_cluster", "CLUSTERER")


		create_clusterer_dag_job_file( cluster_args, get_FINAL_COMMON_ARGS_LINE(), DS_clusterer_outfile_list, \
						       combine_DS_final_outfile, combine_DS_final_condor_file, final_clusterer_num_pose_kept )

		fid_dag.write('\nJOB %s %s\n' % ( combine_DS_final_job_tag, combine_DS_final_condor_file) )

		for clusterer_job_tag in DS_clusterer_job_tag_list:

			if ( clusterer_job_tag not in all_job_tags): error_exit_with_message("clusterer_job_tag (%s) is not in all_job_tags" %(clusterer_job_tag) )
			if ( clusterer_job_tag not in jobs_done ): fid_dag.write('PARENT %s  CHILD %s\n' % (clusterer_job_tag, combine_DS_final_job_tag) )

	#################################Clusterer of ALL_region_FINAL.out##########################################################

	if (exists(ALL_region_FINAL())):

		if (final_rebuild_bulge_step_only==False):

			if (exists(combine_DS_final_outfile)==False):
				error_exit_with_message("ALL_final_outfile (%s) already exist BUT combine_DS_final_outfile (%s) doesn't exist!" %(ALL_region_FINAL(), combine_DS_final_outfile))

			if (exists(standard_region_FINAL())==False):
				error_exit_with_message("ALL_final_outfile (%s) already exist BUT standard_region_FINAL() (%s)  doesn't exist!" %(ALL_region_FINAL(), standard_region_FINAL()))


		print "WARNING: ALL_final_outfile (%s) already exists!" %(ALL_region_FINAL())

	else:

		print "\nCREATE %s JOB" %(ALL_final_job_tag())

		ALL_final_condor_file = get_condor_submit_file(ALL_final_job_tag()+ "_cluster", "CLUSTERER")


		create_clusterer_dag_job_file( cluster_args, get_FINAL_COMMON_ARGS_LINE(), [standard_region_FINAL(), combine_DS_final_outfile], \
						       ALL_region_FINAL(), ALL_final_condor_file, final_clusterer_num_pose_kept )

		fid_dag.write('\nJOB %s %s\n' % ( ALL_final_job_tag(), ALL_final_condor_file) )

		if (exists(combine_DS_final_outfile)==False):
			fid_dag.write('PARENT %s  CHILD %s\n' % (combine_DS_final_job_tag, ALL_final_job_tag()) )

		if (exists(standard_region_FINAL())==False):
			fid_dag.write('PARENT %s  CHILD %s\n' % (standard_region_FINAL(), ALL_final_job_tag()) )

	################################################################################################################################

	print "INCLUDING COMBINE_DS_STEPS: "
	print "Total number of DS_samplerer_jobs to run = %d" %(total_DS_samplerer_jobs)
	print "Total number of 'other' jobs to run = %d" %(len( all_job_tags )-len(jobs_done) )
	print
	print "jobs_already_done   = ", jobs_done
	print
	print
	print "all_job_tags        = ", all_job_tags
	print


#FINAL REBUILD BULGE
if (is_release_mode()==False):
    # This feature is still in testing.

	final_outfile=standard_region_FINAL()
	final_job_tag=standard_final_job_tag()

	if (combine_DS_regions):
		final_outfile=ALL_region_FINAL()
		final_job_tag=ALL_final_job_tag()

	if (final_rebuild_bulge_step_only): #Oct 22, 2011: Hacky!

		fid_dag.close()
		submit_subprocess("rm rna_build.dag")
		fid_dag = open( "rna_build.dag", 'w' ) #Clear all existing JOB!
		fid_dag.write("DOT dag.dot\n")

		if (exists(final_outfile)==False): error_exit_with_message("final_outfile (%s) doesn't exist!" %(final_outfile) )

	#############################################################################################################

	final_rebuild_bulge_job_tag="FINAL_REBUILD_BULGE"

	job_submitted=setup_final_rebuild_bulge_dag_job_file(fid_dag, final_outfile, final_rebuild_bulge_job_tag, native_pdb, \
																			 							final_common_args_file(), "SWA", "")

	if ((exists(final_outfile)==False) and job_submitted):
		fid_dag.write('PARENT %s CHILD %s\n' %(final_job_tag, final_rebuild_bulge_job_tag) )

	#############################################################################################################
	if (BMRB_chemical_shift_file!=""):

		if (exists(BMRB_chemical_shift_file)==False): error_exit_with_message("BMRB_chemical_shift_file (%s) doesn't exist!" %(BMRB_chemical_shift_file) )

		final_rebuild_bulge_job_tag_with_chem_shift="FINAL_REBUILD_BULGE_WITH_CHEM_SHIFT"

		job_submitted=setup_final_rebuild_bulge_dag_job_file(fid_dag, final_outfile, final_rebuild_bulge_job_tag_with_chem_shift, native_pdb, \
																				 							final_common_args_file(), "SWA", BMRB_chemical_shift_file)

		if ((exists(final_outfile)==False) and job_submitted):
			fid_dag.write('PARENT %s CHILD %s\n' %(final_job_tag, final_rebuild_bulge_job_tag_with_chem_shift) )


print "-------------------------------------------------------------------------------------------"
print "Successfully RAN: %s" %(list_to_string(START_argv))
print "-------------------------------------------------------------------------------------------"

