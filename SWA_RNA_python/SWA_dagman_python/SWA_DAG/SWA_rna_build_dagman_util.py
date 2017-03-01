#!/usr/bin/env python


######################################################################
from SWA_DAG_util import *

######################################################################

from SWA_dagman_python.parser.SWA_parse_options import parse_options, replace_arg_value
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import *
######################################################################

from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from os import system, popen
import copy
from sets import Set
######################################################################



if (use_new_src_code()):
	EXE = get_rosetta_EXE("swa_rna_main")
else:
	EXE = get_rosetta_EXE("rna_swa_test")

DB = get_rosetta_database_folder()

##################WARNING, this is far from a comprehensive list of python script used by the Stepwise Assembly method!#########
SAMPLING_PRE_PROCESS_SCRIPT 									= get_PYEXE("SWA_DAG/SWA_sampling_pre_process.py")
SAMPLING_PRE_PROCESS_COMBINE_LONG_LOOP_SCRIPT	= get_PYEXE("SWA_DAG/SWA_sampling_combine_long_loop_pre_process.py")
SAMPLING_POST_PROCESS_SCRIPT 									= get_PYEXE("SWA_DAG/SWA_sampling_post_process.py")
SAMPLE_VIRT_RIBOSE_REDUCER_SCRIPT 						= get_PYEXE("SWA_DAG/SWA_sample_virt_ribose_reducer.py")

#################################################################################################################################

#############################################################################################################################

def import_step_list_from_file_FOUR_EDGE(pathway_file, num_elements):

	print_title_text("Enter import_step_list_from_file=%s " %(pathway_file))
	first_line=True

	step_list=[] ## A list of dictionary. The member of the dictionary are: 'i', 'j', 'bulge_move', 'normal_move'

	num_elements_actual=0

	infile=open(pathway_file,"r")
	for line in infile:
		if (first_line==True):
			first_line=False
			continue

		line_list=line.split(',')
#		print "line_list:", line_list

		step = {}
		step['i_1']=int(line_list[0])
		step['j_1']=int(line_list[1])
		step['i_2']=int(line_list[2])
		step['j_2']=int(line_list[3])

		bulge_move_string=line_list[4]
		bulge_move=False
#		print "bulge_move_string", bulge_move_string

		if (bulge_move_string=="true"):
			step['bulge_move']=True
		elif (bulge_move_string=="false"):
			step['bulge_move']=False
		else:
			error_exit_with_message("Invalid bulge_move_string: %s" %(bulge_move_string))

		normal_move_string=line_list[5]
		normal_move=False
#		print "normal_move_string", normal_move_string

		if (normal_move_string=="true"):
			step['normal_move']=True
		elif (normal_move_string=="false"):
			step['normal_move']=False
		else:
			error_exit_with_message("Invalid normal_move: %s" %(normal_move))

		step['num_elements']=num_elements
		step['tag']=create_tag_FOUR_EDGE(step)

		step_list.append(step)


	for step in step_list:
		print step

	print_title_text("Exit import_step_list_from_file=%s " %(pathway_file))

	return step_list

#############################################################################################################################

def import_step_list_from_file(pathway_file, num_elements):

	print_title_text("Enter import_step_list_from_file=%s " %(pathway_file))
	first_line=True

	step_list=[] ## A list of dictionary. The member of the dictionary are: 'i', 'j', 'bulge_move', 'normal_move'

	num_elements_actual=0

	infile=open(pathway_file,"r")
	for line in infile:
		if (first_line==True):
			first_line=False
			continue

		line_list=line.split(',')
#		print "line_list:", line_list

		step = {}
		step['i']=int(line_list[0])
		step['j']=int(line_list[1])

		i_no_mod=step['i']
		j_no_mod=step['j']
		if (j_no_mod<i_no_mod): j_no_mod=j_no_mod+num_elements

		num_elements_current=j_no_mod-i_no_mod+1
		if (num_elements_actual<num_elements_current): num_elements_actual=num_elements_current

		bulge_move_string=line_list[2]
		bulge_move=False
#		print "bulge_move_string", bulge_move_string

		if (bulge_move_string=="true"):
			step['bulge_move']=True
		elif (bulge_move_string=="false"):
			step['bulge_move']=False
		else:
			error_exit_with_message("Invalid bulge_move_string: %s" %(bulge_move_string))

		normal_move_string=line_list[3]
		normal_move=False
#		print "normal_move_string", normal_move_string

		if (normal_move_string=="true"):
			step['normal_move']=True
		elif (normal_move_string=="false"):
			step['normal_move']=False
		else:
			error_exit_with_message("Invalid normal_move: %s" %(normal_move))
		step_list.append(step)

	for step in step_list:
		print "i= ", step['i'], "j= ", step['j'], "bulge_move= ", step['bulge_move'], "normal_move= ", step['normal_move']

	print_title_text("Exit import_step_list_from_file=%s " %(pathway_file))

	return (step_list, num_elements_actual)

#############################################################################################################################
def get_step_list(pathway):

	print_title_text("Enter get_step_list pathway=%s " %(pathway))
	pathway=pathway.replace("(", "")
	pathway=pathway.replace(")", "")


	raw_steps_list=string.split(pathway,'-')
	print raw_steps_list

	##############Get num_elements##################
	num_elements=0

	for raw_steps in raw_steps_list:

		for raw_step in string.split(raw_steps,','):
			assert(len(raw_step)==1 or len(raw_step)==3)

			if (len(raw_step)==3):
				assert(raw_step[1]=='B')
				bulge_element_num=int( raw_step[0])
				print "bulge_element_num=%d "  %(bulge_element_num)
				if (num_elements<(bulge_element_num+1)): num_elements=bulge_element_num+1

			element_num=int( raw_step[len(raw_step)-1])
			print "element_num=%d "  %(element_num)
			if (num_elements<(element_num+1)): num_elements=element_num+1

	print "num_elements_check=%d " %num_elements
	##############Get num_elements##################

	step_list=[] ## A list of dictionary. The member of the dictionary are: 'i', 'j', 'bulge_move', 'normal_move'
	old_i=0
	old_j=0
	new_i=0
	new_j=0

	for raw_steps in raw_steps_list:

		count=0
		found_bulge_move=False
		found_normal_move=False
		for raw_step in string.split(raw_steps,','):

			count+=1

			assert(len(raw_step)==1 or len(raw_step)==3)
			element_num=int( raw_step[len(raw_step)-1])

			bulge_move=False
			normal_step=False
			step_size=0

			if (len(raw_step)==3):
				assert(raw_step[1]=='B')
				bulge_move=True
				normal_move=False
				step_size=2
			else:
				bulge_move=False
				normal_move=True
				step_size=1

			if (found_bulge_move==False and bulge_move==True): found_bulge_move=True

			if (found_normal_move==False and normal_move==True): found_normal_move=True

			if ((old_i-step_size) % num_elements == element_num):
				new_i  = element_num

				step = {}
				step['i']=new_i
				step['j']=old_j
				step['bulge_move']=bulge_move
				step['normal_move']=normal_move
				step_list.append(step)

			if ((old_j+step_size) % num_elements == element_num):
				new_j  =element_num

				step = {}
				step['i']=old_i
				step['j']=new_j
				step['bulge_move']=bulge_move
				step['normal_move']=normal_move
				step_list.append(step)

		if (count==2):
			step = {}
			step['i']=new_i
			step['j']=new_j
			step['bulge_move']=found_bulge_move
			step['normal_move']=found_normal_move
			step_list.append(step)

		old_i=new_i
		old_j=new_j

	assert( old_i % num_elements == old_j % num_elements )

	for step in step_list:
		print "i= ", step['i'], "j= ", step['j'], "bulge_move= ", step['bulge_move'], "normal_move= ", step['normal_move']

	print_title_text("Exit get_step_list pathway=%s " %(pathway))
	return (step_list, num_elements)
#############################################################################################################################


def satisfy_build_from_helix_to_helix_mode(i, j, num_elements, user_input_pdb_info_list, allow_bulge_right_next_to_input_helix ):

	if (len(user_input_pdb_info_list)>2):
		error_exit_with_message("satisfy_build_from_helix_to_helix_mode called but len(user_input_pdb_info_list)=(%d)>2" %(len(user_input_pdb_info_list) ) )

	modeled_elements_prev = get_modeled_elements(i, j, num_elements)

	num_input_element=0;
	for modeled_element in modeled_elements_prev:

		for pdb_info in user_input_pdb_info_list:
			if (pdb_info["element"]==modeled_element):
				num_input_element += 1
				break

	if (num_input_element>1):
		if ( (j + 1 ) % num_elements == (i) % num_elements ): return True #Close chain

		#this one is for building bulge with single nucleotide step.
		if ( ( (j + 2 ) % num_elements == (i) % num_elements ) and allow_bulge_right_next_to_input_helix ): return True #One step before close chain...

		return False


	return True


##############################################################################################################################

def get_modeled_elements( i, j , num_elements):
    modeled_elements = []
    count = i
    modeled_elements.append( count )
    while not( count == j ):
        count += 1
        count = count % num_elements  #num_elements is defined here???
        modeled_elements.append( count )
    return modeled_elements #This is the list of ALL elements between and including element i and j (in a modolus sense)


########################This function return the residue (belonging to element i) at the boundary between element i and element i+1 #######################################
########################This function assumes num_elements=>3 and that each element occupy a single continuous domain (allow modulo)#######################################
def get_boundary_res( i, assigned_element ):  #assigned_element = []  #map from seq_num to assigned_element
	sample_res = -1
	seq_length = len( assigned_element )
	num_elements = max( assigned_element ) + 1
	target_element = i % num_elements
	target_element_next = (i+1) % num_elements
	for k in range( seq_length  ):
		if (assigned_element[k] == target_element) and (assigned_element[ (k+1) % seq_length ] == target_element_next ):
			sample_res = k+1 #plus+1 since k starts from 0.
			break

	return (sample_res)

#############################################################################################################
def get_job_tag(i, j):
	job_tag = 'REGION_%d_%d' % (i, j)
	return job_tag

#############################################################################################################
def create_tag_FOUR_EDGE(step):

	num_elements=step['num_elements']

	L=len(get_modeled_elements( step['i_1'], step['j_1'] , num_elements) + get_modeled_elements( step['i_2'], step['j_2'] , num_elements))

	tag='L%d_%d_%d_%d_%d' % (L, step['i_1'],step['j_1'], step['i_2'],step['j_2'])

	return tag
#############################################################################################################

def get_job_tag_FOUR_EDGE(step):
		job_tag = 'REGION_' + create_tag_FOUR_EDGE(step)
		return job_tag


#############################################################################################################

def get_region_tags(user_input_pdb_info_list):
		region_tags = []
		for pdb_num in range(len(user_input_pdb_info_list)):
			element = user_input_pdb_info_list[pdb_num]['element']
			region_tags.append(get_job_tag(element, element))
		return region_tags

#############################################################################################################

def satisfy_optimize_long_loop_mode(optimize_long_loop_mode, OLLM_chain_closure_only, i, j, i_prev, j_prev, num_elements):

		if (optimize_long_loop_mode == False): return False

		if ( i <=  j  ): return False;

		# Ok no bulge move when combining the two sides for now!
		if ( ( i_prev % num_elements) - ( i% num_elements) ) % num_elements == 2:
			return False

		if ( ( j % num_elements) - ( j_prev % num_elements) ) % num_elements == 2:
			return False

		if i == 0:
			return False
		if j == 0:
			return False


		if OLLM_chain_closure_only:
			# Only include steps right at the chain-closure.
			if i % num_elements != (j + 1) % num_elements:
				return False

			if i_prev % num_elements != (j_prev + 2) % num_elements:
				return False

		return True

#############################################################################################################
def screen_dinucleotide_at_single_element_chain_closure( dinucleotide_at_single_element_cc, element_definition, user_input_pdb_info_list, i, j, moving_element, num_elements):

	Is_valid_value_dinucleotide_at_single_element_cc(dinucleotide_at_single_element_cc)

	if ( ( i% num_elements) != ( (j+1) % num_elements )  ): return True

	if ( len(element_definition[moving_element])!=1 ): return True

	if (dinucleotide_at_single_element_cc=="all"):
		return True
	elif (dinucleotide_at_single_element_cc=="pure_append_prepend"):
		if (i==0 or j==0): return True
	elif (dinucleotide_at_single_element_cc=="none"):
		return False
	else:
		error_exit_with_message("Invalid dinucleotide_at_single_element_cc value=%s" %(dinucleotide_at_single_element_cc) )

	return False

################################################################################################################################################################

def Is_valid_attachment(i, j, i_prev, j_prev , num_elements, cutpoints_open, all_job_tags, pathway_bulge_res, element_definition, user_input_pdb_info_list, \
														optimize_long_loop_mode, OLLM_chain_closure_only, allow_bulge_right_next_to_input_helix, \
														allow_normal_move_at_pathway_bulge_res, allow_bulge_move_at_non_pathway_bulge_res, \
														dinucleotide_at_single_element_cc, build_from_helix_to_helix_mode, check_attach_verbose):


		if (check_attach_verbose): print "---Check is valid attachment: i=%d, i_prev=%d, j=%d, j_prev=%d--- " %( i, i_prev, j, j_prev )

		attach_side=""

		bulge_element=9999
		moving_element=9999
		attach_location=9999


		num_possible_attach_side=0
		Is_dinucleotide_step=False

		if (i!=i_prev):
			attach_side="prepend"
			num_possible_attach_side+=1
			moving_element=i

			if ( ( (i_prev-i) % num_elements ) == 1 ):
				Is_dinucleotide_step=False
			elif ( ( (i_prev-i) % num_elements ) == 2):
				bulge_element=(i + 1) % num_elements
				Is_dinucleotide_step=True
			else:
				error_exit_with_message( "(i_prev-i) % num_elements ) != 1 or 2, (i_prev-i) % num_elements )=%d" %((i_prev-i) % num_elements ) )


		if (j!=j_prev):
			attach_side="append"
			num_possible_attach_side+=1
			moving_element=j

			if ( ( (j-j_prev) % num_elements ) == 1 ):
				Is_dinucleotide_step=False
			elif ( ( (j-j_prev) % num_elements ) == 2):
				bulge_element=(j - 1) % num_elements
				Is_dinucleotide_step=True
			else:
				error_exit_with_message( "(i_prev-i) % num_elements ) != 1 or 2, (i_prev-i) % num_elements )=%d" %((i_prev-i) % num_elements ) )


		if (num_possible_attach_side!=1): error_exit_with_message( "num_possible_attach_side=(%d)!=1" %(num_possible_attach_side) )


		prev_job_tag = get_job_tag(i_prev,j_prev)

		can_attach=True

		if ( ( prev_job_tag not in all_job_tags) and satisfy_optimize_long_loop_mode(optimize_long_loop_mode, OLLM_chain_closure_only,  i, j, i_prev, j_prev, num_elements)==False ):
			can_attach=False

		if (attach_side=="prepend"):
			if ( element_definition[moving_element][-1] in cutpoints_open): can_attach=False #if cutpont_open at 5' edge then cannot prepend
		elif (attach_side=="append"):
			if ( (element_definition[moving_element][0]-1 ) in cutpoints_open): can_attach=False #if cutpont_open at 3' edge then cannot append
		else:
			error_exit_with_message( "Invalid attach_side=%s" %(attach_side) )


		if (build_from_helix_to_helix_mode):
			if ( satisfy_build_from_helix_to_helix_mode(i, j, num_elements, user_input_pdb_info_list, allow_bulge_right_next_to_input_helix )==False  ):  can_attach=False


		if (Is_dinucleotide_step):

			#if (allow_bulge_right_next_to_input_helix==False): #Mod out on May 31, 2011. Floating base + lenghth(moving_element)>1 is not implemented yet!
			if ( len(element_definition[moving_element])!=1 ): can_attach=False


			if ( len(element_definition[bulge_element])!=1 ): can_attach=False

			if ( 	len(pathway_bulge_res)!=0 and (allow_bulge_move_at_non_pathway_bulge_res==False) ):

				if (attach_side=="prepend"):
					if (element_definition[bulge_element][-1] not in pathway_bulge_res): can_attach=False
				else:
					if (element_definition[bulge_element][0]  not in pathway_bulge_res): can_attach=False

			if (screen_dinucleotide_at_single_element_chain_closure(dinucleotide_at_single_element_cc, element_definition, user_input_pdb_info_list, i, j, moving_element, num_elements)==False): can_attach=False

		else:

			if ( ( ( i% num_elements) != ( (j+1) % num_elements ) )  and allow_normal_move_at_pathway_bulge_res==False):

				if (attach_side=="prepend"):
					if ( element_definition[moving_element][-1] in pathway_bulge_res): can_attach=False
				else:
					if ( element_definition[moving_element][0]  in pathway_bulge_res): can_attach=False


		if (check_attach_verbose):
			print "Is_valid_attachment=%s, moving_element=%d, bulge_element=%d, attach_side=%s " %(can_attach, moving_element, bulge_element, attach_side)

		return can_attach



#############################################################################################################
def get_silent_file_pose_info(i_prev, j_prev, num_elements, element_definition):

		modeled_elements_prev = get_modeled_elements( i_prev, j_prev , num_elements)
		if ( len( modeled_elements_prev ) == 1 ):
			error_exit_with_message("len( modeled_elements_prev ) == 1!  i_prev = %d, j_prev = %d " % (i_prev, j_prev) )

		infile = 'region_%d_%d_sample.cluster.out' % (i_prev,j_prev)
		tag = 'S_$(Process)' #name of pose in the silent_file.

		modeled_res_prev = []
		for m in modeled_elements_prev:
			for k in element_definition[ m ]: modeled_res_prev.append( k ) #element_definition[ m ] does not have to be ordered
		modeled_res_prev.sort()

		pdb_info={}
		pdb_info['tag'] = tag
		pdb_info['residue_list']=modeled_res_prev
		pdb_info['silent_file']=infile
		pdb_info['i_prev']=i_prev
		pdb_info['j_prev']=j_prev

		return pdb_info
#############################################################################################################

def get_modeled_pose_info_list_FOUR_EDGE(prev_step, num_elements, base_pdb_tag, act_input_pdb, element_definition):

	#	print "Enter get_modeled_pose_info_list_FOUR_EDGE function"
		prev_job_tag=get_job_tag_FOUR_EDGE(prev_step)

		modeled_elements_prev = get_modeled_elements( prev_step['i_1'], prev_step['j_1'] , num_elements) + get_modeled_elements( prev_step['i_2'], prev_step['j_2'] , num_elements)

		modeled_res_prev = []
		for m in modeled_elements_prev:
			for k in element_definition[ m ]: modeled_res_prev.append( k ) #element_definition[ m ] does not have to be ordered
		modeled_res_prev.sort()

		pdb_info={}
		pdb_info['residue_list']=modeled_res_prev


		if prev_job_tag in base_pdb_tag:

			pdb_info['tag'] = act_input_pdb
			pdb_info['silent_file'] = ''

		else:

			if ( len( modeled_elements_prev ) == 1 ):  error_exit_with_message("len( modeled_elements_prev ) == 1 for prev_step " + get_job_tag_FOUR_EDGE(rev_step) )

			infile = 'region_%s_sample.cluster.out' % ( prev_step['tag'].lower())
			tag = 'S_$(Process)' #name of pose in the silent_file.


			pdb_info['tag'] = tag
			pdb_info['silent_file']=infile

		#print pdb_info
		return [pdb_info]
#############################################################################################################
def get_modeled_pose_info_list_combine_long_loop(i, j, i_prev, j_prev, num_elements, user_input_pdb_info_list, element_definition, optimize_long_loop_mode, OLLM_chain_closure_only):

		input_pose_info_list = []

		input_pdb_tags=get_region_tags(user_input_pdb_info_list)

		if (satisfy_optimize_long_loop_mode(optimize_long_loop_mode, OLLM_chain_closure_only, i, j, i_prev, j_prev, num_elements)):

			print 'satisfy_optimize_long_loop_mode(i=%d,j=%d, i_prev=%d, j_prev=%d, num_elements=%d)==true' %(i,j,i_prev,j_prev, num_elements)

			##note that the order below is important, since the c++ assume that the first input silent_file is APPEND and the second input silent_file is PREPEND

			if get_job_tag(0, j_prev) not in input_pdb_tags:	 #APPEND
				pdb_info_1=get_silent_file_pose_info(0, j_prev, num_elements, element_definition)
				input_pose_info_list.append(pdb_info_1)

			if get_job_tag(i_prev, 0) not in input_pdb_tags: #PREPEND
				pdb_info_2=get_silent_file_pose_info(i_prev, 0, num_elements, element_definition)
				input_pose_info_list.append(pdb_info_2)




		return input_pose_info_list


#############################################################################################################

def get_modeled_pose_info(i, j, i_prev, j_prev, num_elements, build_from_scratch_tags, user_input_pdb_info_list, element_definition):


		pose_info={}
		prev_job_tag=get_job_tag(i_prev, j_prev)
		input_pdb_tags=get_region_tags(user_input_pdb_info_list)

		if prev_job_tag in input_pdb_tags:

			pdb_num = input_pdb_tags.index( prev_job_tag )
			pose_info = user_input_pdb_info_list[pdb_num]

		elif ( prev_job_tag in build_from_scratch_tags):

			fake_input_element=0
			if ( i == i_prev ):
				fake_input_element = i
			elif ( j == j_prev):
				fake_input_element = j
			else:
				error_exit_with_message("(i != i_prev) and  (j != j_prev)")

			if ( len( element_definition[ fake_input_element ] ) != 1 ):
				print 'fake_input_element= %d, len( element_definition[ fake_input_element ] ), %d ' %( fake_input_element, len( element_definition[ fake_input_element ] ) )
				error_exit_with_message("len( element_definition[ fake_input_element ] ) != 1")

			pose_info['tag'] = "build_from_scratch"
			pose_info['residue_list']=[ element_definition[ fake_input_element ][0] ]
			pose_info['silent_file']=""
			pose_info['i_prev']=i_prev
			pose_info['j_prev']=j_prev

		else:
			pose_info=get_silent_file_pose_info(i_prev, j_prev, num_elements, element_definition)

		return pose_info


#############################################################################################################

def get_user_input_pdb_info_list(input_pdbs, input_silent_file_list, input_res_full, sequence):

		for i in range( len(input_res_full) ): assert( input_res_full[i] not in input_res_full[:i] )

		if (input_silent_file_list!=[""]):
			if (len(input_pdbs)!=len(input_silent_file_list)):
				error_exit_with_message("input_silent_file_list!=[""] but len(input_pdbs)!=len(input_silent_file_list) !!")

		user_input_pdb_info_list = []

		count = 0
		input_res  = [] #Order of elements in input_res matches the order of elements in input pdb
		if ( len(input_pdbs[0])> 0 ):
			pdb_num=0
			for input_pdb in input_pdbs:

				if ( exists( input_pdb )==False ): error_exit_with_message("exists( %s )==False" %(input_pdb) )

				#input_pdb_sequence = popen( "pdb2fasta.py %s "  %(input_pdb) ).readlines()[1][:-1]
				input_pdb_sequence = popen( "get_sequence.py %s "  %(input_pdb) ).readlines()[0][:-1]
				input_seq_length = len( input_pdb_sequence )

				starting_sequence = ""
				input_res_particular = []
				for i in range( count, count+input_seq_length ):
					starting_sequence += sequence[ input_res_full[ i ] - 1  ]
					input_res_particular.append( input_res_full[ i ] )

				if (starting_sequence != input_pdb_sequence) :
					error_exit_with_message('start_sequence =%s does not equal input_pdb_sequence=%s ' % (starting_sequence, input_pdb_sequence))

				count += input_seq_length

				input_pdb_info = {}


				if (input_silent_file_list!=[""] and input_silent_file_list[pdb_num]!="NONE"):

					if (exists(input_silent_file_list[pdb_num])==False):
						error_exit_with_message("input_silent_file_list[%d] %s doesn't exist!" %(pdb_num, input_silent_file))
					input_pdb_info['tag'] = 'S_$(Process)'
					input_pdb_info['silent_file'] = input_silent_file_list[pdb_num]

				else:

					input_pdb_info['tag'] = input_pdb
					input_pdb_info['silent_file'] = ''

				input_pdb_info['residue_list'] = input_res_particular
				input_pdb_info['element'] = -1

				user_input_pdb_info_list.append(input_pdb_info)
				pdb_num+=1

			print "user_input_pdb_info_list: "
			output_pdb_info_list(user_input_pdb_info_list)

		return user_input_pdb_info_list

#############################################################################################################
def output_pdb_info_list(pdb_info_list):
		for pdb_info in pdb_info_list:
#			print pdb_info
#			print 'tag: %s  silent_file: %s  element: %d residue_list: ' % (pdb_info['tag'], pdb_info['silent_file'], pdb_info['element']), pdb_info['residue_list']
			print 'tag: %s  silent_file: %s ' % (pdb_info['tag'], pdb_info['silent_file'])

#############################################################################################################

def satisfy_duplex_frame_slip_util(i,j, max_duplex_frame_slip, duplex_frame_BP):

	high_element=0
	low_element=0


	if (j<i): #build inward
		high_element=i
		low_element=j
	else: #build_outward
		high_element=j
		low_element=i

	diff_low_element=low_element-duplex_frame_BP[0]
	diff_high_element=high_element-duplex_frame_BP[1]

	if ( abs(diff_low_element+diff_high_element)> max_duplex_frame_slip):
		return False
	else:
		return True

#############################################################################################################


def satisfy_duplex_frame_slip(i,j, num_elements, max_duplex_frame_slip, duplex_frame_BP):

	if (i==j): return True

	if (len(duplex_frame_BP)==0):
		if (max_duplex_frame_slip!=0): error_exit_with_message("len(duplex_frame_BP)==0 but max_duplex_frame_slip!=0" )
		return True

	if (len(duplex_frame_BP)!=2): error_exit_with_message("User input in duplex_frame_BP but len(duplex_frame_BP)!=2" )

	if (max_duplex_frame_slip<=0): error_exit_with_message("max_duplex_frame_slip<=0" )

	if (duplex_frame_BP[0]>=duplex_frame_BP[1]):
		error_exit_with_message("duplex_frame_BP[0]>=duplex_frame_BP[1], duplex_frame_BP=%s " %(list_to_string(duplex_frame_BP) ) )

	#max_duplex_frame_slip= parse_options(argv, "max_duplex_frame_slip", 0)
	#duplex_frame_BP= parse_options(argv, "max_duplex_frame_slip", [ 0 ]) # these are elements

	if (satisfy_duplex_frame_slip_util(i,j, max_duplex_frame_slip, duplex_frame_BP)): return True

	if (i==0):
		if (satisfy_duplex_frame_slip_util(num_elements,j, max_duplex_frame_slip, duplex_frame_BP)): return True

	if (j==0):
		if (satisfy_duplex_frame_slip_util(i, num_elements, max_duplex_frame_slip, duplex_frame_BP)): return True

	return False



################################################################################################################

def	satisfy_enforce_path_base_pair(enforce_path_base_pair_list, i, i_prev, j, j_prev, allow_bulge, num_elements):

	element_list=get_modeled_elements(i, j, num_elements)

	prev_element_list=get_modeled_elements(i_prev, j_prev, num_elements)

	possible_moving_element=0
	if (i!=i_prev):
		moving_element=i
		possible_moving_element+=1

	if (j!=j_prev):
		moving_element=j
		possible_moving_element+=1

	if (possible_moving_element!=1): error_exit_with_message("possible_moving_element!=1")

	allow_offset=1
	if (allow_bulge): allow_offset=2

	for base_pair in enforce_path_base_pair_list:


		if (base_pair[0] not in element_list and base_pair[1] not in element_list): continue #Haven't modeled either base!

		if (base_pair[0] in element_list): #check that base_pair[0] is not skipped
			if (base_pair[0] not in 	prev_element_list and (moving_element!=base_pair[0]) ):
				print "prev_element_list=", prev_element_list
				print "(i=%d, j=%d, i_prev=%d, j_prev=%d BLAH1.1" %(i, j, i_prev, j_prev)
				return False

		if (base_pair[1] in element_list): #check that base_pair[1] is not skipped
			if (base_pair[1] not in 	prev_element_list and (moving_element!=base_pair[1]) ):
				print "(i=%d, j=%d, i_prev=%d, j_prev=%d BLAH1.2" %(i, j, i_prev, j_prev)
				return False

		if (base_pair[0] in element_list and base_pair[1] in element_list): continue #Already modeled both base!

		fail_enforce_path_base_pair=True

		#OK..now only one of the two BP elements is in modeled_element_list.
		if (base_pair[0] in element_list):

			if (base_pair[1] in element_list): error_exit_with_message("base_pair[1] in element_list")

			if (Is_subset_list(get_modeled_elements(i, base_pair[0], num_elements) , get_modeled_elements(base_pair[1], base_pair[0], num_elements) ) ):

				if ( j==base_pair[0] ):
					if (j!=j_prev): #build j
						#build base_pair[0] (j) only if base_pair[1] (i) can be built at the next step.
						if ( len(element_list)+allow_offset >= len(get_modeled_elements(base_pair[1], base_pair[0], num_elements))	):  fail_enforce_path_base_pair=False
					if (i!=i_prev): #build i
						if ( i==base_pair[1] ): fail_enforce_path_base_pair=False


			if (Is_subset_list(get_modeled_elements(base_pair[0], j, num_elements) , get_modeled_elements(base_pair[0], base_pair[1], num_elements) ) ):

				if ( i==base_pair[0] ):
					if (i!=i_prev): #build i
						#build base_pair[0] (i) only if base_pair[1] (j) can be built at the next step.
						if ( len(element_list)+allow_offset >= len(get_modeled_elements(base_pair[0], base_pair[1], num_elements))	): fail_enforce_path_base_pair=False
					if (j!=j_prev): #build j
						if ( j==base_pair[1] ): fail_enforce_path_base_pair=False

		##########################################################################################
		if (base_pair[1] in element_list):

			if (base_pair[0] in element_list): error_exit_with_message("base_pair[0] in element_list")

			if (Is_subset_list(get_modeled_elements(base_pair[1], j, num_elements), get_modeled_elements(base_pair[1], base_pair[0], num_elements) ) ):

				if ( i==base_pair[1] ):
					if (i!=i_prev): #build i
						#build base_pair[1] (i) only if base_pair[0] (j) can be built at the next step.
						if ( len(element_list)+allow_offset >= len(get_modeled_elements(base_pair[1], base_pair[0], num_elements))	): fail_enforce_path_base_pair=False
					if (j!=j_prev): #build j
						if ( j==base_pair[0] ): fail_enforce_path_base_pair=False

			if (Is_subset_list(get_modeled_elements(i, base_pair[1], num_elements), get_modeled_elements(base_pair[0], base_pair[1], num_elements) ) ):

				if ( j==base_pair[1] ):
					if (j!=j_prev): #build j
						#build base_pair[1] (j) only if base_pair[0] (i) can be built at the next step.
						if ( len(element_list)+allow_offset >= len(get_modeled_elements(base_pair[0], base_pair[1], num_elements))	): fail_enforce_path_base_pair=False
					if (i!=i_prev): #build i
						if ( i==base_pair[0] ): fail_enforce_path_base_pair=False

		if (fail_enforce_path_base_pair==True): return False

	return True
#############################################################################################################
'''
def	satisfy_enforce_path_base_pair_OLD(enforce_path_base_pair_list, i, i_prev, j, j_prev, allow_bulge, num_elements):

	modeled_element_list=get_modeled_elements(i, j, num_elements)

	fail_enforce_path_base_pair=False
	for base_pair in enforce_path_base_pair_list:




		i_act=i
		i_prev_act=i_prev
		j_act=j

		if (i==0): i_act=num_elements
		if (i_prev==0): i_prev_act=num_elements

		if ( i==0 and (( i% num_elements) == ( (j+1) % num_elements )) ): #hacky, Nov 16, 2010
			i_act=0

		if ( j==0 and (( i% num_elements) == ( (j+1) % num_elements )) ): #hacky, Nov 16, 2010
			j_act=num_elements

		i_BP=0
		j_BP=0

		if (i_act>j_act):
			i_BP=base_pair[1]
			j_BP=base_pair[0]
		elif (j_act>i_act):
			i_BP=base_pair[0]
			j_BP=base_pair[1]
		else:
			error_exit_with_message("i_act(%d)==j_act(%d)!" %(i_act, j_act) )

		if (j_act<j_BP and i_act<i_BP): fail_enforce_path_base_pair=True
		if (i_act>i_BP and j_act>j_BP): fail_enforce_path_base_pair=True

		if (j_act!=j_prev and j_act==j_BP): #build j and j is a base pair.
			if (i_act!=i_BP and i_act!=i_BP+1 and (i_act!=i_BP+2 and allow_bulge) ):  fail_enforce_path_base_pair=True
			#Build j only if i has been built or can be build at the next step

		if (j_act!=j_prev and i==i_BP): #build j and i is a base pair
			if ( j_act<j_BP ):  fail_enforce_path_base_pair=True  #Build j only if j_BP has been built or if j==j_BP
			if ( j_act>j_BP and j_prev<j_BP ): fail_enforce_path_base_pair=True #check that j_BP is not skipped


		if (i!=i_prev and i==i_BP): #build i and i is a base pair
			if ( j_act!=j_BP and j_act!=j_BP-1 and (j_act!=j_BP-2 and allow_bulge) ):  fail_enforce_path_base_pair=True
			#Build i only if j has been built or can be build at the next step

		if (i!=i_prev and j_act==j_BP): #build i and j is a base pair.
			if (i_act>i_BP ): fail_enforce_path_base_pair=True #base_pair[1] is being built or has already been built
			if ( i_act<i_BP and i_prev_act>i_BP ): fail_enforce_path_base_pair=True #check that i_BP is not skipped

	if (fail_enforce_path_base_pair==True):
		return False
	else:
		return True
'''

#############################################################################################################
def extract_silent_file_info(pose_info_list):

		silent_files=""
		silent_file_tags=""
		num_silent_files=0

		for pose_info in pose_info_list:
			if (pose_info["silent_file"]!=""):
				num_silent_files+=1
				silent_files+="%s " %(pose_info["silent_file"])
				silent_file_tags+="%s " %(pose_info["tag"])

		return (num_silent_files, silent_files, silent_file_tags)

#############################################################################################################

def get_num_silent_files(pose_info_list):

		(num_silent_files, silent_files, silent_file_tags)=extract_silent_file_info(pose_info_list)

		return num_silent_files



#############################################################################################################

def get_filterer_job_tag( sampler_job_tag):

	filterer_job_tag='%s_FILTERER' % ( sampler_job_tag )

	return filterer_job_tag

#############################################################################################################


def get_filterer_outfile( sampler_job_tag):

	filterer_outfile='FILTER_OUTFILE/%s_output.txt' %( get_filterer_job_tag( sampler_job_tag).lower() )

	return filterer_outfile


#############################################################################################################
def add_input_res_to_common_args(job_specific_common_args, pose_info_list):
		if (len(pose_info_list)>2 ):
			output_pdb_info_list(pose_info_list)
			error_exit_with_message("length of pose_info_list (%d) is greater than 2!" %(len(pose_info_list)) )

		for n in range(len(pose_info_list)):
			pose_info=pose_info_list[n]
			if (n==0):
				job_specific_common_args += ' -input_res '
			else:
				job_specific_common_args += ' -input_res2 '

			job_specific_common_args += list_to_string_with_dashes( pose_info[ "residue_list" ] )


		return job_specific_common_args

#############################################################################################################

def	create_samplerer_dag_job_file(job_specific_common_args, job_specific_sampling_args, job_specific_filterer_args, \
																filterer_mode, pose_info_list, output_foldername, job_tag, prefix):

		#print "enter_create_samplerer_dag_job_file, job_tag=%s" %(job_tag)
		###########################add input_pose_info to sampling and filterer_arg#################

		if (len(pose_info_list)>2 ):
			output_pdb_info_list(pose_info_list)
			error_exit_with_message("length of pose_info_list (%d) is greater than 2!" %(len(pose_info_list)) )


		(num_silent_files, silent_files, silent_file_tags)=extract_silent_file_info(pose_info_list)


		pdb_tags=""

		for pose_info in pose_info_list:
			if (pose_info["silent_file"]==""): pdb_tags+="%s " % (pose_info["tag"])

		if (pdb_tags!=""):     job_specific_sampling_args += ' -s %s ' %(pdb_tags)
		if (silent_files!=""):
			if (num_silent_files==2): #need filterer.
				job_specific_sampling_args += ' -in:file:silent %s ' %(silent_files)
				job_specific_sampling_args += ' -job_queue_ID $(Process) ' #determine -tags from job_queue_ID
				job_specific_filterer_args += ' -in:file:silent %s ' %(silent_files)
			else:
				job_specific_sampling_args += ' -in:file:silent %s -tags %s ' %(silent_files, silent_file_tags)


		################################################################################################

		queue_tag = 'S_$(Process)'  ###umm this specific form is assumed in SWA_sampling_pre_process.py

		if (num_silent_files==0): #only one job in the queue, pre-process script is not called
			output_foldername_act=output_foldername + '/ONLY_JOB'
		else:
			output_foldername_act=output_foldername + '/' + queue_tag

		outfile = output_foldername_act + '/' + prefix + 'sample.out'
		job_specific_sampling_args += '-out:file:silent %s ' % outfile

		mapper_outfiles_sampling=' ' + outfile


		##################Submit the filterer job########################
		if ( num_silent_files==2):

			filterer_job_tag= get_filterer_job_tag( job_tag )
			filterer_outfile= get_filterer_outfile( job_tag )

			if (exists( dirname(filterer_outfile) )==False): submit_subprocess( 'mkdir -p %s' %(dirname(filterer_outfile)) )

			job_specific_filterer_args += ' -filter_output_filename %s ' %(filterer_outfile)
			job_specific_sampling_args += ' -filter_output_filename %s ' %(filterer_outfile) #will need to open this file and read in the tag..

			if (filterer_mode=="combine_DS_regions"):

				job_specific_filterer_args += ' -combine_helical_silent_file true  '
				job_specific_sampling_args += ' -combine_helical_silent_file true  '

			elif (filterer_mode=="combine_long_loop_allow_previous_clash"):

				#### filtering for previous clashes is causing problems in native
				#### rmsd_screen runs -- specifically, 5P_j55a_group_I_intron
				#### CHANGE turn off filter_for_previous_clash -- will revert if causes problems
				#### - Caleb April 10, 2015      
				job_specific_filterer_args += ' -filter_for_previous_clash false  '
				job_specific_sampling_args += ' -combine_long_loop_mode true '

			elif (filterer_mode=="combine_long_loop_clash"):

				job_specific_filterer_args += ' -filter_for_previous_clash true  '
				job_specific_sampling_args += ' -combine_long_loop_mode true '

			elif (filterer_mode=="combine_long_loop_contact"):

				job_specific_filterer_args += ' -filter_for_previous_contact true '
				job_specific_sampling_args += ' -combine_long_loop_mode true '

			else:
				error_exit_with_message("num_silent_files==2 but invalid filterer_mode=%s " %(filterer_mode) )

			filterer_condor_submit_file = get_condor_submit_file( filterer_job_tag , "FILTERER")

			make_dag_job_submit_file( filterer_condor_submit_file, EXE, job_specific_filterer_args + ' ' + job_specific_common_args, '', filterer_outfile, '', '')

		##################Submit the filterer job########################

		condor_submit_file = get_condor_submit_file( job_tag, "SAMPLER" )

		#print "creating dag_job_submit_file=%s " %(job_tag)

		make_dag_job_submit_file( condor_submit_file, EXE, job_specific_sampling_args + ' ' + job_specific_common_args, '', mapper_outfiles_sampling, '', get_sampler_post_process_filter_outfile(output_foldername) )


####################################################################################################

def submit_sampling_pre_and_post_process(fid_dag, job_tag, pose_info_list, output_foldername, min_post_process_filtered_struct):

		filterer_outfile= get_filterer_outfile( job_tag )
		condor_submit_file = get_condor_submit_file( job_tag , "SAMPLER")

		(num_silent_files, silent_files, silent_file_tags)=extract_silent_file_info(pose_info_list)

		##################PRE_PROCESS#################

		if (num_silent_files==2): fid_dag.write('SCRIPT PRE %s   %s %s %s  \n' % (job_tag, SAMPLING_PRE_PROCESS_COMBINE_LONG_LOOP_SCRIPT , filterer_outfile ,  condor_submit_file ) )

		if (num_silent_files==1): fid_dag.write('SCRIPT PRE %s   %s %s %s \n' % (job_tag, SAMPLING_PRE_PROCESS_SCRIPT  , silent_files , condor_submit_file ) )

		##################POST_PROCESS#################


		fid_dag.write('SCRIPT POST %s %s -indir_prefix %s -condor_submit_file %s -min_filtered_struct %d\n' %(job_tag, SAMPLING_POST_PROCESS_SCRIPT, output_foldername, \
																																																		 condor_submit_file, min_post_process_filtered_struct))

####################################################################################################


def update_dag_dependency( fid_dag, job_tag, pose_info_list, filterer_parent_job_tag_list, sampler_parent_job_tag_list, all_job_tags, jobs_done):

	#(parent_job_tag not in jobs_done) --> #check that parent_job_tag is not accomplished in a prior BIOX run or an user input element

	if ( get_num_silent_files(pose_info_list)==2):

		filterer_job_tag=get_filterer_job_tag( job_tag )

		filterer_condor_submit_file = get_condor_submit_file( filterer_job_tag, "FILTERER" )

		for filterer_parent_job_tag in filterer_parent_job_tag_list:

			if (filterer_parent_job_tag not in all_job_tags): error_exit_with_message("filterer_parent_job_tag (%s) is not in all_job_tags!" %(filterer_parent_job_tag) )

			if (filterer_parent_job_tag not in jobs_done): fid_dag.write('PARENT %s  CHILD %s\n' %(filterer_parent_job_tag, filterer_job_tag) )

		fid_dag.write('\nJOB %s %s\n' % (filterer_job_tag, filterer_condor_submit_file) )
		fid_dag.write('PARENT %s  CHILD %s\n' % (filterer_job_tag, job_tag) )

	#######################################################################################
	fid_dag.write('\nJOB %s %s\n' % (job_tag, get_condor_submit_file( job_tag, "SAMPLER" ) ) )

	for sampler_parent_job_tag in sampler_parent_job_tag_list:

		if (sampler_parent_job_tag not in all_job_tags): error_exit_with_message("sampler_parent_job_tag (%s) is not in all_job_tags!" %(sampler_parent_job_tag) )
		if (sampler_parent_job_tag not in jobs_done): fid_dag.write('PARENT %s  CHILD %s\n' % (sampler_parent_job_tag, job_tag) )



#######################################################################################

def 	create_clusterer_dag_job_file( cluster_args, job_specific_common_args, input_silent_file_list, clusterer_outfile, condor_submit_cluster_file, clusterer_num_pose_kept	):

		if (exists( clusterer_outfile ) ): error_exit_with_message("clusterer_outfile %s already exist!" %(clusterer_outfile) )

		for n in range(len(input_silent_file_list)): #check that silent_file in the list is not a empty string!
			if (basename(input_silent_file_list[n])==""): error_exit_with_message("basename(input_silent_file_list[%d])==\"\"!" %(n))

		specific_cluster_args = cluster_args + job_specific_common_args + ' -in:file:silent %s -out:file:silent %s  -clusterer_num_pose_kept %d ' % (string.join( input_silent_file_list ), clusterer_outfile , clusterer_num_pose_kept)

		make_dag_job_submit_file( condor_submit_cluster_file, EXE, specific_cluster_args, '', clusterer_outfile, '', '')

#######################################################################################

#Reorder the seq_num from prepend attachment point to append attachment point (from the perspective of the element itself)
#Note this doesn't correspond to normal low to high integer sorting espeically for outer element
#THIS IMPLICITLY ASSUME THAT THERE ARE TWO ATTACHMENT EDGES.
def reorder_element_definition(element_definition, cutpoints_open, num_elements, total_res):

	if (len(element_definition)!=num_elements): error_exit_with_message("len(element_definition)!=num_elements" )


	for element_id in range( num_elements ):

		seq_list=copy.deepcopy(element_definition[ element_id ])

		seq_list.sort() #normal sorting from low to high (IN_PLACE)

		seq_break_list=[]

		for n in range( len( seq_list) ):
			if (n==0): continue

			curr_seq_num=	seq_list[n]
			prev_seq_num=	seq_list[n-1]

			if ( (prev_seq_num+1)!=curr_seq_num ): seq_break_list.append(prev_seq_num)

		if (len(seq_break_list)>1):
				print "seq_list= ", seq_list ," seq_break_list= ", seq_break_list
				error_exit_with_message("len(seq_break_list)>1" )


		#if (len(seq_list)>1):
		#	found_CP=False

		#	for CP_open in cutpoints_open:
		#		if (CP_open in seq_list): found_CP=True

		#	if (found_CP==False): Is_outer_element=True

		Is_outer_element=False

		if (len(seq_break_list)==1): Is_outer_element=True

		print "element_id=%d, Is_outer_element=%s" %(element_id, Is_outer_element), "seq_list=" , seq_list

		if (Is_outer_element):

			if (seq_list[0]!=1): #for now, the only outer_element is outer DS helix containing the first and last res.
				error_exit_with_message("seq_list[0]!=1, seq_list[0]=%d" %(seq_list[0]) )

			if (total_res!=seq_list[-1]): #for now, the only outer_element is outer DS helix containing the first and last res.
				error_exit_with_message("total_res!=seq_list[-1], total_res=%d, seq_list[-1]=%d" %(total_res, seq_list[-1]) )

			reordered_seq_list=[]


			for seq_num in seq_list:
				if (seq_num>seq_break_list[0]): reordered_seq_list.append(seq_num)

			for seq_num in seq_list:
				if (seq_num<=seq_break_list[0]): reordered_seq_list.append(seq_num)

			seq_list=reordered_seq_list

		else:
			if (len(seq_break_list)!=0):
				error_exit_with_message("Is_outer_element==False but len(seq_break_list)!=0" )

		element_definition[ element_id ]=seq_list

	return element_definition


###########################################################################
def get_combine_DS_REGION_FOLDER():
	return "COMBINE_DS_REGIONS/"

#June 13, 2011
def get_create_input_ds_combine_region_job_tag(region_ID):

	job_tag='CREATE_INPUT_DS_COMBINE_REGION_%d_%d' % (region_ID[0],region_ID[1])

	return job_tag

########################################################################
def create_input_silent_file_for_DS_combine(region_ID, pose_info, cluster_args, all_job_tags, jobs_done, create_input_silent_for_DS_combine_job_tags, fid_dag, NSTRUCT, global_sample_res_plus_bound):

	job_tag=get_create_input_ds_combine_region_job_tag(region_ID)

	combine_DS_folder=get_combine_DS_REGION_FOLDER()

	input_DS_silent_file=combine_DS_folder + "/INPUT_SILENT_FILES/input_DS_combine_region_%d_%d.out" % (region_ID[0], region_ID[1])

	if (job_tag not in all_job_tags): all_job_tags.append(job_tag)

	create_dag=True

	###NEED THIS? Since DAG_continuous.py does also check if the job is already done DURING RUNTIME!####
	###ACTUALLY NEED THIS, since DAG_continuous.py only deals with incompleted job where the reducer done file doesn't exist!
	if (exists(input_DS_silent_file)): #Already done this job in a prior SWA_rna_build_dagman.py run.
		jobs_done.append(job_tag)
		create_dag=False

	if (job_tag in create_input_silent_for_DS_combine_job_tags): #Already submitted this job in the current SWA_rna_build_dagman.py run
		create_dag=False

	if (create_dag):

		create_input_silent_for_DS_combine_job_tags.append(job_tag)

		common_args=safe_readlines( 'COMMON_ARGS/common_args_region_%d_%d_sample.cluster.out' % (region_ID[0],region_ID[1]) )[0]

		##FEB 10, 2012: I did some testing for the GAGU duplex with removing these next two lines, but decided to revert back!
		common_args=replace_arg_value(common_args, "alignment_res", [list_to_string(global_sample_res_plus_bound,'-')[1:]],  allow_missing=True)
		common_args=replace_arg_value(common_args, "filter_user_alignment_res", "false",  allow_missing=True)

		MOD_cluster_args=cluster_args
		MOD_cluster_args=replace_arg_value(MOD_cluster_args, "suite_cluster_radius", 2.0, allow_missing=False)
		MOD_cluster_args=replace_arg_value(MOD_cluster_args, "loop_cluster_radius", 1.5,  allow_missing=False)
		MOD_cluster_args=replace_arg_value(MOD_cluster_args, "clusterer_quick_alignment", "false",  allow_missing=True)

		input_DS_condor_submit = get_condor_submit_file( job_tag , "DS_COMBINE_CLUSTERER")

		if (basename(pose_info['silent_file'])==""): error_exit_with_message("basename(pose_info['silent_file'])==\"\"")

		create_clusterer_dag_job_file( MOD_cluster_args, common_args, [pose_info['silent_file']], input_DS_silent_file, input_DS_condor_submit , NSTRUCT )

		fid_dag.write('\nJOB %s %s\n' % ( job_tag, input_DS_condor_submit) )

		standard_sampling_clusterer_job_tag = get_job_tag(region_ID[0],region_ID[1])

		if ( standard_sampling_clusterer_job_tag not in all_job_tags): error_exit_with_message("standard_sampling_clusterer_job_tag (%s) is not in all_job_tags" %(standard_sampling_clusterer_job_tag) )
		if ( standard_sampling_clusterer_job_tag not in jobs_done ): fid_dag.write('PARENT %s  CHILD %s\n' % (standard_sampling_clusterer_job_tag, job_tag) )

	return input_DS_silent_file

##################################################################################
def get_virt_ribose_sampler_folder():
	return "VIRTUAL_RIBOSE_SAMPLER/"

##################################################################################
def get_virt_ribose_sampler_job_tag(pose_info, standard_sampler_job_tag):

	job_tag='VIRT_RIBOSE_SAMPLER_REGION_%d_%d_FOR_%s' % (pose_info['i_prev'],pose_info['j_prev'], standard_sampler_job_tag)

	return job_tag

##################################################################################
def create_sampled_virt_ribose_silent_file(pose_info, standard_sampler_job_tag, prev_clusterer_job_tag, all_job_tags, jobs_done, sample_virtual_ribose_list, native_pdb, fid_dag, num_silent_files):
	#TODO (sripakpa): explain why have to sample virt ribose for each seperate build step that uses a particular outfile.

	verbose=True

	if (pose_info['i_prev']==pose_info['j_prev']): return (pose_info['silent_file'], True) #No virtual ribose in user inputted region!

	sample_virt_ribose_folder=get_virt_ribose_sampler_folder()

	job_tag=get_virt_ribose_sampler_job_tag(pose_info, standard_sampler_job_tag)

	if (job_tag in all_job_tags): error_exit_with_message("job_tag (%s) already in all_job_tags!" %(job_tag) )
	all_job_tags.append(job_tag)

	common_args=safe_readlines( 'COMMON_ARGS/common_args_region_%d_%d_sample.cluster.out' % (pose_info['i_prev'],pose_info['j_prev']) )[0]

	job_condor_submit = get_condor_submit_file( job_tag.replace("VIRT_RIBOSE_SAMPLER_", "") , "VIRT_RIBOSE_SAMPLER")

	####Check if job is done####
	log_foldername = get_DAG_log_foldername(job_condor_submit)

	reducer_job_info = get_job_info(log_foldername, 'reducer')

	queue_tag = 'S_$(Process)'  ###umm this specific form is assumed in SWA_sampling_pre_process.py

	output_foldername="%s/%s/" %(sample_virt_ribose_folder, job_tag.replace("VIRT_RIBOSE_SAMPLER_", ""))

	reducer_outfile="%s/reducer_sample_virt_ribose_region_%d_%d.out" % ( output_foldername, pose_info['i_prev'], pose_info['j_prev'])

	DAG_already_done=False

	#OK DAG is done! Important, cannot submit DAG that already done, since DAG_continuous.py will raise error if reducer_done_signal_filename already exist!
	if ( exists(reducer_job_info['done_signal_filename']) ):
		if (exists(reducer_outfile)==False): error_exit_with_message("reducer_outfile (%s) doesn't exist!" %(reducer_outfile) )
		DAG_already_done=True

	if (DAG_already_done==False):

		if (basename(pose_info['silent_file'])==""): error_exit_with_message("basename(pose_info['silent_file']")

		if (verbose): print "standard_sampler_job_tag= %s, sample_virtual_sugar_list %s " %(standard_sampler_job_tag, list_to_string(sample_virtual_ribose_list) )

		#####################Simple preprocess script (same as preprocess of standard sampler)######################################################
		# this has to go after job declaration in condor DAGMAN file... moved below
		#fid_dag.write('SCRIPT PRE %s   %s %s %s\n' % (job_tag, SAMPLING_PRE_PROCESS_SCRIPT  , pose_info['silent_file'] , job_condor_submit) )

		########################################################################################################################

		mapper_output_foldername=output_foldername + queue_tag

		mapper_outfile = mapper_output_foldername + '/sample_virt_ribose_region_%d_%d.out' %(pose_info['i_prev'], pose_info['j_prev'])

		job_args = ' -algorithm rna_sample_virtual_sugar -database %s ' %(DB)
		job_args += common_args
		job_args += ' -in:file:silent %s -tags %s -out:file:silent %s ' %(pose_info['silent_file'], pose_info["tag"], mapper_outfile)
		job_args += ' -sample_virtual_sugar_list %s ' %(list_to_string(sample_virtual_ribose_list))

		if ( native_pdb!="" ):
			if (exists( native_pdb )==False): error_exit_with_message("exists( native_pdb (%s) )==False" %(native_pdb))
			job_args += '-native %s ' %(native_pdb)

		make_dag_job_submit_file( job_condor_submit, EXE, job_args, '', mapper_outfile, '', reducer_outfile)

		fid_dag.write('\nJOB %s %s\n' % ( job_tag, job_condor_submit) )

		# moved from above. need PRE processing script defined after JOB declaration
		fid_dag.write('SCRIPT PRE %s   %s %s %s \n' % (job_tag, SAMPLING_PRE_PROCESS_SCRIPT  , pose_info['silent_file'] , job_condor_submit ) )

		if ( prev_clusterer_job_tag not in all_job_tags): error_exit_with_message("prev_clusterer_job_tag (%s) is not in all_job_tags" %(standard_sampler_job_tag) )
		if ( prev_clusterer_job_tag not in jobs_done ): fid_dag.write('PARENT %s  CHILD %s\n' % (prev_clusterer_job_tag, job_tag) )

		###################Post process, combine the structures into a single silent_file#######################################
		post_process_command= SAMPLE_VIRT_RIBOSE_REDUCER_SCRIPT
		post_process_command+=" -data_foldername %s " %(output_foldername)
		post_process_command+=" -reducer_outfile %s " %(basename(reducer_outfile) )
		post_process_command+=" -START_silent_file %s " %(pose_info['silent_file'])

		need_combine_silent_file_filterer=(num_silent_files==2) #Kinda hacky, is there a better way to determine this?

		#Need this for the combine_silent_file_filterer to work for -filterer_undercount_ribose_rotamers mode!
		if (need_combine_silent_file_filterer):
			post_process_command+=" -keep_old_score True "
			post_process_command+=" -keep_old_tag True "

		fid_dag.write('SCRIPT POST %s %s\n' % (job_tag, post_process_command ) )


	return ( reducer_outfile, DAG_already_done )

#############################################################################################################

def standard_region_FINAL():

	return "region_FINAL.out"

#############################################################################################################
def standard_final_job_tag():

	return "REGION_FINAL"

#############################################################################################################
def ALL_region_FINAL():

	return "ALL_region_FINAL.out"

#############################################################################################################
def ALL_final_job_tag():

	return "ALL_REGION_FINAL"

#############################################################################################################

def final_common_args_file():

	return 'COMMON_ARGS/common_args_region_FINAL.out'

#############################################################################################################

def get_FINAL_COMMON_ARGS_LINE():

	if (exists(final_common_args_file())==False): error_exit_with_message("FINAL_COMMON_ARGS (%s) doesn't exist!" %(final_common_args_file()))

	line_list=safe_readlines(final_common_args_file())

	if ( len(line_list) !=1 ): error_exit_with_message("# len(line_list)=(%s)!=1" %( len(line_list) ) )

	if ( line_list[0][-1]=='\n'): error_exit_with_message("line_list[0][-1]=='\n'")

	return line_list[0]

#############################################################################################################
