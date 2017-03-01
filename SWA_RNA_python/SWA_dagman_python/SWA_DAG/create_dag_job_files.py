#!/usr/bin/env python

######################################################################
from SWA_dagman_python.utility.DAGMAN_util import create_generic_README_SUB
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

if(is_release_mode()==False): #Not yet RELEASED!
	from SWA_dagman_python.utility.DAGMAN_util import setup_chemical_shift_args
######################################################################
from sys import argv,exit
import sys
import traceback
from glob import glob
import string
from os.path import basename, dirname, exists, expanduser
from time import sleep
from sets import Set
######################################################################

start_argv=copy.deepcopy(argv)

#create_dag_job_files.py

#create_dag_job_files.py -native_pdb native_CA_BP.pdb -native_rmsd_screen False -nstruct 1000 -num_slave_nodes 250 -cutpoint_open  5 -protonated_H1_adenosine_list 8 > LOG_create_dag_job_files.out 2> LOG_create_dag_job_files.err

#create_dag_job_files.py -native_pdb native_CA_BP.pdb -native_rmsd_screen False -nstruct 1000 -num_slave_nodes 250 -cutpoint_open  5  > LOG_create_dag_job_files.out 2> LOG_create_dag_job_files.err

#create_dag_job_files.py -fasta 2r8s_fasta_with_GUA_loop -s  remove_AA_Trimmed_2r8s_RNA.pdb -sample_res 11 12 13 -native_virtual_res 12 -single_stranded_loop_mode True

#create_FARFAR_job_files.py  -s  mutate_C123-G143_remove_AA_2r8s_RNA.pdb -sample_res 123 124 125 -native_virtual_res 124 -single_stranded_loop_mode True > create_SWA_job_files_LOG.txt

#create_FARFAR_job_files.py -fasta 2r8s_fasta_with_GUA_loop -s  remove_AA_Trimmed_2r8s_RNA.pdb -sample_res 11 12 13 -native_virtual_res 12 -single_stranded_loop_mode True

#create_dag_job_files.py -native_rmsd_screen False  -native_pdb FINAL_April_2010_C7_2_tetraloop_receptor_model.pdb -s mutated_no_loop_2BP_upper_helix_Trimmed_2r8s_RNA.pdb -nstruct 1000  -sample_res 11 12 13 -native_virtual_res 12 -num_slave_nodes 500 -single_stranded_loop_mode True > create_dag.out 2> create_dag.err

#-------------------------------------------------------------KINK TURN-----------------------------------------------------------------------------------------------#

#SHORT
#create_dag_job_files.py -native_rmsd_screen False  -native_pdb kink_turn_3iqn.pdb -nstruct 1000 -native_virtual_res 12 -cutpoint_open  7  -num_slave_nodes 500  > create_dag.out 2> create_dag.err

#LONG
#create_dag_job_files.py -native_rmsd_screen False  -native_pdb long_kink_turn_3iqn.pdb -nstruct 1000 -native_virtual_res 14 -cutpoint_open  8  -num_slave_nodes 500  > create_dag.out 2> create_dag.err

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#


################################################################################
single_stranded_loop_mode=parse_options(argv, "single_stranded_loop_mode", "False") #long_loop_mode

fasta_file=parse_options( argv, "fasta", "" )

if(is_release_mode() and single_stranded_loop_mode):

	rna_torsion_potential_folder=parse_options(argv, "rna_torsion_potential_folder", "ps_03242010/" )
	force_field_file=parse_options(argv, "force_field_file", 'stepwise/rna/rna_loop_hires_04092010.wts') #Context Independent Geometric Solvation but no intra_base_phosphate terms.
	sample_virt_ribose_in_sep_DAG=parse_options(argv, "sample_virt_ribose_in_sep_DAG", "True" )
	apply_VDW_rep_delete_matching_res=parse_options(argv, "apply_VDW_rep_delete_matching_res", "True")
	clusterer_optimize_memory_usage=parse_options(argv, "clusterer_optimize_memory_usage", "true" )

	if(fasta_file==""): error_exit_with_message("User need to specified fasta option!")

else:

	rna_torsion_potential_folder=parse_options(argv, "rna_torsion_potential_folder", "" )
	force_field_file=parse_options(argv, "force_field_file", 'stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts')
	sample_virt_ribose_in_sep_DAG=parse_options(argv, "sample_virt_ribose_in_sep_DAG", "True" )
	apply_VDW_rep_delete_matching_res=parse_options(argv, "apply_VDW_rep_delete_matching_res", "False")
	clusterer_optimize_memory_usage=parse_options(argv, "clusterer_optimize_memory_usage", "false" )


################################################################################

native_pdb= parse_options( argv, "native_pdb", "" )
cutpoint_open=parse_options( argv, "cutpoint_open", 0)
num_slave_nodes=parse_options( argv, "num_slave_nodes", 500)
native_rmsd_screen=parse_options( argv, "native_rmsd_screen", "False")
rmsd_screen=parse_options( argv, "rmsd_screen", 0.0 )
clusterer_num_pose_kept=parse_options(argv, "nstruct", 1000)
native_virtual_res=parse_options(argv, "native_virtual_res", [-1])
force_bulge_res=parse_options(argv, "force_bulge_res", [-1])
sample_res_list=parse_segment_string_list( parse_options(argv, "sample_res", [""]) )

clusterer_keep_pose_in_memory=parse_options(argv, "clusterer_keep_pose_in_memory", "true" )

old_SWA_idealize_helix=parse_options(argv, "old_SWA_idealize_helix", "False" ) #For Backward compatibility after switch to new_idealized helix on Feb 03, 2012

allow_combine_DS_regions=parse_options( argv, "allow_combine_DS_regions", "False")

no_bulge = parse_options( argv, "no_bulge", "False")

OLLM_allow_previous_clash=parse_options(argv, "OLLM_allow_previous_clash", "False")

tether_jump = parse_options( argv, "tether_jump", "True")

VDW_rep_optimize_memory_usage = parse_options( argv, "VDW_rep_optimize_memory_usage", "False")

#################################################################################
VDW_rep_screen_info_list=parse_options(argv, "VDW_rep_screen_info", [""])

if(VDW_rep_screen_info_list==[""]): VDW_rep_screen_info_list=[]

#################################################################################
start_elements=parse_options(argv, "s", [""])
if(start_elements==[""]): start_elements=[]
print "start_elements= ", start_elements
for pdb in start_elements:
	if(exists(pdb)==False): error_exit_with_message("start_pdb (%s) doesn't exist!!" %(pdb))
#################################################################################


extra_anti_chi_rotamer=parse_options(argv, "extra_anti_chi_rotamer", "false")
extra_syn_chi_rotamer=parse_options(argv, "extra_syn_chi_rotamer", "false")

force_syn_chi_res_list=parse_options(argv, "force_syn_chi_res_list", [-1] )
force_north_ribose_list=parse_options(argv, "force_north_ribose_list", [-1] )
force_south_ribose_list=parse_options(argv, "force_south_ribose_list", [-1] )
protonated_H1_adenosine_list=parse_options(argv, "protonated_H1_adenosine_list", [-1] )
enforce_path_base_pairs=parse_options(argv, "enforce_path_base_pairs", [""] ) #these are elements pairs (e.g 3-8 5-6)

#################################################################################
if(is_release_mode()==False): #Not yet RELEASED!
	ignore_chemical_shift=parse_options( argv, "ignore_chemical_shift", "False")

	BMRB_chemical_shift_file=parse_options(argv, "BMRB_chemical_shift_file", "" )

	if(ignore_chemical_shift):
		BMRB_chemical_shift_file=""
		BLAHBLAHBLAH=parse_options( argv, "chemical_shift_H5_prime_mode", "" )
	else:
		setup_chemical_shift_args(argv)
#################################################################################

for seq_num in force_north_ribose_list:
	if(seq_num in force_south_ribose_list): error_exit_with_message("seq_num %s is in both force_north_ribose_list and force_south_ribose_list!" %(seq_num) )

invert_lower_element= parse_options( argv, "invert_lower_element", "False" ) #Nov 23, 2010  this is to ensure prefect alignment
invert_upper_element= parse_options( argv, "invert_upper_element", "False" ) #Nov 23, 2010  this is to ensure prefect alignment

if(invert_lower_element==True):
	print "WARNING, user specified that invert_lower_element==True!"

if(invert_upper_element==True):
	print "WARNING, user specified that invert_upper_element==True!"
#################################################################################

#Between Feb 16, 2011 and March 28, 2011: dinucleotide_at_single_element_cc="pure_append_prepend" for both single_stranded_loop_mode and denovo motif case.
dinucleotide_at_single_element_cc="none"

#Before Feb 16, this used to "none"
if(single_stranded_loop_mode): dinucleotide_at_single_element_cc="pure_append_prepend"

print "dinucleotide_at_single_element_cc= %s " %(dinucleotide_at_single_element_cc)

Is_valid_value_dinucleotide_at_single_element_cc(dinucleotide_at_single_element_cc)

#################################################################################

denovo_modeling=True
if(len(sample_res_list)!=0): denovo_modeling=False  #only modeling part of the RNA..

###################extra FARFAR options that have no meaning in SWA#####
BLAH_BLAH_BLAH=parse_options(argv, "randomize_ribose_puckers", "False", Verbose=False)
BLAH_BLAH_BLAH=parse_options(argv, "fragment_match_type", "", Verbose=False)
BLAH_BLAH_BLAH=parse_options(argv, "allow_bulge_mode", "", Verbose=False)
BLAH_BLAH_BLAH=parse_options(argv, "torsion_database", "", Verbose=False)
BLAH_BLAH_BLAH=parse_options(argv, "excise_segment_torsion_DB", [""], Verbose=False)
BLAH_BLAH_BLAH=parse_options(argv, "allow_consecutive_bulges", "False", Verbose=False)
######################################

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )


verbose=True

if(fasta_file!=""):

	if(exists(fasta_file)==False): error_exit_with_message("User specified fasta_file (%s) doesn't exist!" %(fasta_file))

elif(native_pdb!=""):

	if(exists(native_pdb)==False):  error_exit_with_message("User specified native_pdb (%s) doesn't exist!!" %(native_pdb))

	fasta_file="fasta"

	submit_subprocess('SWA_pdb2fasta.py %s > %s' %(native_pdb, fasta_file))

else:
	error_exit_with_message("User did not specified either fasta_file or native_pdb option!")


#####create the fasta file##########

sequence = open( fasta_file ).readlines()[1][:-1]
if(verbose): print "sequence= %s " %(sequence)

total_res=len(sequence)
all_res=range(1,total_res+1)


####Feb 03, 2012: crate SETUP_LOG folder############################################################
if(exists("SETUP_LOG/")): error_exit_with_message("SETUP_LOG/ already exist!")

submit_subprocess("mkdir SETUP_LOG/")

####create README_SETUP.py ###########################################################
input_res=[]
fixed_res=[]
jump_point_pair_list=[]

num_strands=1
if(cutpoint_open!=0): num_strands=2

#################################################################################################
if(denovo_modeling):

	if(len(start_elements)!=0): error_exit_with_message("denovo_modeling==true but user passed in start_elements(%s) (-s)" %(list_to_string(start_elements) ) )
	if(len(sample_res_list)!=0): error_exit_with_message("denovo_modeling==true but user passed in sample_res_list(%s)" %(list_to_string(sample_res_list) ) )


	###########Feb 03, 2012#######
	####create idealized starting elements and VDW_rep_screener########################################
	if(len(VDW_rep_screen_info_list)!=0): error_exit_with_message("User pass in VDW_rep_screen_info for a denovo modeling case")

	if(exists('VDW_rep_screen_info.txt')==True): error_exit_with_message("VDW_rep_screen_info.txt already exist!")

	if(old_SWA_idealize_helix==False):

		lower_helix_filename="lower_helix.pdb"
		upper_helix_filename="upper_helix.pdb"

		create_lower_helix   ="create_idealize_helix_general.py -sequence %s -verbose %s -VDW_screener_type LOWER " %(sequence, verbose)
		create_lower_helix  +="-strand_1_seq_num %s -strand_2_seq_num %s " %(list_to_string([1,2]), list_to_string([total_res-1, total_res]))
		create_lower_helix  +="-helix_filename %s > SETUP_LOG/LOG_create_lower_helix.txt " %(lower_helix_filename)
		submit_subprocess(create_lower_helix)

		if(num_strands==2):
			create_upper_helix ="create_idealize_helix_general.py -sequence %s -verbose %s -VDW_screener_type UPPER " %(sequence, verbose)
			create_upper_helix+="-strand_1_seq_num %s -strand_2_seq_num %s " %(list_to_string([cutpoint_open-1,cutpoint_open]), list_to_string([cutpoint_open+1,cutpoint_open+2]))
			create_upper_helix+="-helix_filename %s > SETUP_LOG/LOG_create_upper_helix.txt " %(upper_helix_filename)
			submit_subprocess(create_upper_helix)

	else:
		######Deprecated on Feb 03, 2012######
		lower_helix_filename="lower_elements.pdb"
		upper_helix_filename="upper_elements.pdb"

		submit_subprocess("create_idealize_helix.py -fasta %s -verbose %s -invert_helix %s " %(fasta_file, verbose, invert_lower_element) )

		if(num_strands==2):
			submit_subprocess("create_idealize_helix.py -fasta %s -cutpoint_open %d -verbose %s -invert_helix %s " %(fasta_file, cutpoint_open, verbose, invert_upper_element) )
		######Deprecated on Feb 03, 2012######


	if(exists('VDW_rep_screen_info.txt')==False): error_exit_with_message("VDW_rep_screen_info.txt doesn't exist!")

	VDW_rep_screen_info_lines=open( "VDW_rep_screen_info.txt" ).readlines()

	if(len(VDW_rep_screen_info_lines)!=num_strands):
		error_exit_with_message( "len(VDW_rep_screen_info_lines)=(%s)!=(%s)=num_strands" %(len(VDW_rep_screen_info_lines), num_strands) )

	for n in range(len(VDW_rep_screen_info_lines)):
		print "VDW_rep_screen_info_lines[%d]=%s" %(n, VDW_rep_screen_info_lines[n]),
		line_list=VDW_rep_screen_info_lines[n].split()

		if(line_list[0]!="-VDW_rep_screen_info"): error_exit_with_message("line_list[0]!=\"-VDW_rep_screen_info\"")

		VDW_rep_screen_info_list.extend(line_list[1:])

	submit_subprocess("rm VDW_rep_screen_info.txt")

	##################################

	start_elements=[lower_helix_filename]
	if(num_strands==2): start_elements.append(upper_helix_filename)

	input_res=[1,2,total_res-1, total_res]
	if(num_strands==2): input_res.extend([cutpoint_open-1,cutpoint_open,cutpoint_open+1,cutpoint_open+2])

	fixed_res=[1, total_res]
	if(num_strands==2): fixed_res.extend([cutpoint_open,cutpoint_open+1])

	jump_point_pair_list=['%d-%d' %(fixed_res[0], fixed_res[1])]
	if(num_strands==2): jump_point_pair_list.append('%d-%d' %(fixed_res[2],fixed_res[3]) )

	sample_res_list=list( Set(all_res)-Set(input_res) )

else: #long loop modeling and etc.
	if(len(start_elements)==0): error_exit_with_message("denovo_modeling==false but user did not pass in start_elements (-s)")
	if(len(sample_res_list)==0): error_exit_with_message("denovo_modeling==False but len(sample_res_list)==0!")

	input_res=list( Set(all_res)-Set(sample_res_list) )
	input_res.sort()

	fixed_res=list( Set(all_res)-Set(sample_res_list) ) #for long loop mode
	fixed_res.sort()

	jump_point_pair_list=['%d-%d' %(all_res[0], all_res[-1])]  #for long loop mode

#################################################################################################
if(len(VDW_rep_screen_info_list)>0): check_valid_VDW_rep_screen_info_list(VDW_rep_screen_info_list)

sample_res_list.sort() #make sure the list is sorted.


native_alignment_res=[]
for seq_num in all_res: #basically the sample_res and its intermediate neighboring res.

	if( (seq_num in sample_res_list) and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)
	if( ((seq_num+1) in sample_res_list) and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)
	if( ((seq_num-1) in sample_res_list) and (seq_num not in native_alignment_res) ): native_alignment_res.append(seq_num)

native_alignment_res.sort()

if(single_stranded_loop_mode): native_alignment_res=[] #Jan 7, 2010...don't use native_alignment in long loop mode since non-loop region is static and align exactly with native_pdb.

alignment_res=jump_point_pair_list[:]

rmsd_res=list(Set(all_res)-Set(input_res) ) #Change on Oct 17, calculate rmsd only over building res.
rmsd_res.sort()

sample_segment_list=get_segment_string_list(sample_res_list)


print "     sample_res_list= ", sample_res_list
print " sample_segment_list= ", sample_segment_list
print "            rmsd_res= ", rmsd_res
print "native_alignment_res= ", native_alignment_res


README_SETUP = open( "README_SETUP.py", 'w')

README_SETUP.write( '#!/usr/bin/env python\n' )
README_SETUP.write( 'from os import system\n' )
README_SETUP.write( 'import string\n\n' )

README_SETUP.write( "command= '%s '\n\n" %(get_PYEXE("~/src/rosetta/tools/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_build_dagman.py")) )
README_SETUP.write( "command+= '-s %s -fasta %s '\n\n" %(list_to_string(start_elements), fasta_file) )

if(native_pdb!=""): README_SETUP.write( "command+= '-native %s '\n\n" %(native_pdb) )

README_SETUP.write( "command+= '-clusterer_num_pose_kept %s '\n\n" %(clusterer_num_pose_kept) )

if(cutpoint_open!=0):
	README_SETUP.write( "command+= '-cutpoint_open %d '\n\n" %(cutpoint_open ) )

README_SETUP.write( "command+= '-force_field_file %s '\n\n" %(force_field_file))

README_SETUP.write( "command+= '-sampler_extra_anti_chi_rotamer %s '\n\n" %(extra_anti_chi_rotamer) ) #change to false on Oct 21, 2010...speed up code.
README_SETUP.write( "command+= '-sampler_extra_syn_chi_rotamer %s '\n\n"  %(extra_syn_chi_rotamer)  ) #change to false on Oct 21, 2010...speed up code.

README_SETUP.write( "command+= '-clusterer_quick_alignment true '\n\n" ) #NOT default option of SWA_parse_rosetta_options!
README_SETUP.write( "command+= '-clusterer_optimize_memory_usage %s '\n\n" %(clusterer_optimize_memory_usage) )

#Dec 8, 2010...if keep just 1000 pose, that two_stage mode actually take up more space and is actually much slower
README_SETUP.write( "command+= '-clusterer_keep_pose_in_memory %s '\n\n" %(clusterer_keep_pose_in_memory) )

README_SETUP.write( "command+= '-allow_combine_DS_regions %s '\n\n" %(allow_combine_DS_regions) ) #Feb 07, 2012

#############################################

README_SETUP.write( "command+= '-dinucleotide_at_single_element_cc %s '\n\n" %(dinucleotide_at_single_element_cc) ) #implement this on Feb 16, 2011

#	sampling_argv+= " -PBP_clustering_at_chain_closure true "  Should turn this on? Oct 23, 2010

if(len(native_virtual_res) > 0):
	README_SETUP.write( "command+= '-native_virtual_res %s '\n\n" %(list_to_string(native_virtual_res))  )

if( ( native_rmsd_screen ) or ( rmsd_screen > 0.0 )):
	if( native_rmsd_screen ): README_SETUP.write( "command+= '-native_rmsd_screen true '\n\n" )
	if( rmsd_screen > 0.0  ): README_SETUP.write( "command+= '-rmsd_screen %d '\n\n" % (rmsd_screen) )

	README_SETUP.write( "command+= '-sampler_num_pose_kept 40 '\n\n" )

	if(len(native_virtual_res) > 0):
		README_SETUP.write( "command+= '-pathway_bulge_res %s '\n\n" %(list_to_string(native_virtual_res))  )
		README_SETUP.write( "command+= '-allow_bulge True '\n\n")
	else:
		README_SETUP.write( "command+= '-allow_bulge False '\n\n")
else: #Full exploration of conformational space
	if ( no_bulge ):
		README_SETUP.write( "command+= '-allow_bulge False '\n\n")
	else:
		README_SETUP.write( "command+= '-allow_bulge True '\n\n")

	if(len(force_bulge_res)>0):
		assert( not no_bulge )
		README_SETUP.write( "command+= '-pathway_bulge_res %s '\n\n" %(list_to_string(force_bulge_res))  )


README_SETUP.write( "command+= '-input_res %s '\n\n" %(list_to_string(input_res) ) )
README_SETUP.write( "command+= '-rmsd_res %s '\n\n" %(list_to_string(rmsd_res) ) )

if(len(native_alignment_res)!=0):
	README_SETUP.write( "command+= '-native_alignment_res %s '\n\n" %(list_to_string(native_alignment_res) ) )

README_SETUP.write( "command+= '-jump_point_pair_list %s '\n\n" %(list_to_string(jump_point_pair_list) ) )
README_SETUP.write( "command+= '-fixed_res %s '\n\n" %(list_to_string(fixed_res) ) )
README_SETUP.write( "command+= '-alignment_res_list %s '\n\n" %(list_to_string(alignment_res) ) )

if(len(VDW_rep_screen_info_list) > 0):
	README_SETUP.write( "command+= '-VDW_rep_screen_info %s '\n\n" %(list_to_string(VDW_rep_screen_info_list) ) )

	if(len(sample_segment_list)==0): error_exit_with_message("len(sample_segment_list)==0")

	if(apply_VDW_rep_delete_matching_res):
		README_SETUP.write("command+= '-VDW_rep_delete_matching_res %s '\n\n" %(list_to_string(sample_segment_list) ) )
	else:
		README_SETUP.write("command+= '-VDW_rep_delete_matching_res false '\n\n" )


if(single_stranded_loop_mode):
	README_SETUP.write( "command+= '-optimize_long_loop_mode True '\n\n" )
	README_SETUP.write( "command+= '-OLLM_chain_closure_only True '\n\n" )
	if OLLM_allow_previous_clash:
		README_SETUP.write( "command+= '-OLLM_allow_previous_clash True '\n\n" )
else:
	README_SETUP.write( "command+= '-analytic_etable_evaluation False '\n\n" ) 


if(len(force_syn_chi_res_list)>0): README_SETUP.write( "command+= '-force_syn_chi_res_list %s '\n\n" %(list_to_string(force_syn_chi_res_list) ) )

if(len(force_south_ribose_list)>0): README_SETUP.write( "command+= '-force_south_ribose_list %s '\n\n" %(list_to_string(force_south_ribose_list) ) )

if(len(force_north_ribose_list)>0): README_SETUP.write( "command+= '-force_north_ribose_list %s '\n\n" %(list_to_string(force_north_ribose_list) ) )

if(enforce_path_base_pairs!=[""]): README_SETUP.write( "command+= '-enforce_path_base_pairs %s '\n\n" %(list_to_string(enforce_path_base_pairs) ) )

if(len(protonated_H1_adenosine_list)>0): README_SETUP.write( "command+= '-protonated_H1_adenosine_list %s '\n\n" %(list_to_string(protonated_H1_adenosine_list) ) )

if(rna_torsion_potential_folder!=""): README_SETUP.write( "command+= '-rna_torsion_potential_folder %s '\n\n" %(rna_torsion_potential_folder) )

if(is_release_mode()==False): #Not yet RELEASED!
	if(BMRB_chemical_shift_file!=""): README_SETUP.write( "command+= '-BMRB_chemical_shift_file %s '\n\n" %(BMRB_chemical_shift_file) )

if ( not no_bulge ): README_SETUP.write( "command+= '-sample_virt_ribose_in_sep_DAG %s '\n\n" %(sample_virt_ribose_in_sep_DAG) )

if ( no_bulge ): README_SETUP.write( "command+= '-floating_base False -allow_bulge_at_chainbreak false ' \n\n" )

if ( not tether_jump ):
	README_SETUP.write( "command+= '-tether_jump false ' \n\n" )

if ( VDW_rep_optimize_memory_usage ):
	README_SETUP.write( "command+= '-VDW_rep_optimize_memory_usage true ' \n\n" )

README_SETUP.write( "command+= '>LOG_SWA_rna_build_dagman.out '\n\n" )

README_SETUP.write( "system(command) \n\n" )

README_SETUP.close()

#disallow options
#command+= '-pathway_file specific_pathway.txt '
#command+= '-allow_build_from_scratch True '
#command+= '-allow_bulge_right_next_to_input_helix False ####still need to implement this for floating base mode..'


####create README_SUB ############################################################
create_generic_README_SUB(num_slave_nodes)

####################################################################################
print "----------------------------------------------------------------------------------------------------------------------------"
print "Successfully RAN %s" %( list_to_string(start_argv) )
print "----------------------------------------------------------------------------------------------------------------------------"
