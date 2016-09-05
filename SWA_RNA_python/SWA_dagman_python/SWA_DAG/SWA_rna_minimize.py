#!/usr/bin/env python

######################################################################

from SWA_rna_minimize_util import *
######################################################################

from SWA_dagman_python.utility.DAGMAN_util import create_generic_README_SUB
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option, get_option_name_args_safe, ensure_no_duplicate_options

######################################################################

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy

######################################################################


#plot '<grep "syn_chi_torsion"  output.txt ' u 2:4
#syn_chi_torsion= -25.4147 syn_chi_offset= 24.5853

#plot '<grep "ratio   5"   swa_rna_minimize_local.out  ' u 11:12, x
#plot '<grep "ratio   5"   swa_rna_minimize_local.out  ' u 11:($12-$11), 0
#plot '<grep "ratio   5"   swa_rna_minimize_local.out  ' u 11:(100*($12-$11)/$12), 100, -100, 0 ls 3

#plot '<grep "ratio   5"   rna_minimize_output.txt  ' u 11:12, x
#plot '<grep "ratio   5"   srna_minimize_output.txt  ' u 11:($12-$11), 0
#plot '<grep "ratio   5"   rna_minimize_output.txt  ' u 11:(100*($12-$11)/$12), 100, -100, 0 ls 3




#core.optimization: ratio inc rsd typ atm prsd ptyp patm natm   numeric  analytic     ratio       f11       f00       f22  vars[ii]
#numeric--> column 11
#analytic--> column 12

#########################################################
###BIOX2###:
#~/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_minimize.py -silent_file WITH_SHIFT_STATS_Jan_23_SWA.out   -force_field_file rna/rna_infinitesimal_CS_03022011_w_intra_terms.wts   -chemical_shift_file renumbered_short_bmr5705.bmrb    -chemical_shift_res_list DEFAULT   -common_args common_args_region_FINAL.out -native_pdb short_2koc_RNA_A.pdb -num_slave_nodes 500    -LOCATION BIOX2  -remove_virtual_res_variant_during_minimize true  > SWA_rna_minimize_LOG.txt 

#-deriv_check true 

#~/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_minimize.py -silent_file WITH_SHIFT_STATS_March_04_SWA.out   -force_field_file rna/rna_4X_CS_03022011_w_intra_terms.wts   -chemical_shift_file renumbered_MANUAL_SRP_28MER_EXP_shift.bmrb   -chemical_shift_res_list DEFAULT  -common_args common_args_region_FINAL.out -native_pdb 1lnt_fix_non_natural_base.pdb -num_slave_nodes 500    -LOCATION BIOX2 -remove_virtual_res_variant_during_minimize true > SWA_rna_minimize_LOG.txt 

###UCAC_tetraloop
#~/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_minimize.py -silent_file WITH_SHIFT_STATS_Nov_10_FARFAR.out  -force_field_file rna/rna_4X_CS_03022011_w_intra_terms.wts    -chemical_shift_file renumbered_AUTO_UCAC_shift.bmrb   -chemical_shift_res_list DEFAULT  -common_args common_args_region_FINAL.out -native_pdb G1C8_UCAC_1s72.pdb -num_slave_nodes 500    -LOCATION BIOX2 -remove_virtual_res_variant_during_minimize true > SWA_rna_minimize_LOG.txt 

#WITH_SHIFT_STATS_Feb_06_SWA.out 
# -chemical_shift_file adjusted_AUTO_UCAC_shift.bmrb
#-chemical_shift_file renumbered_AUTO_UCAC_shift.bmrb  -tag_range 0 1000 
#-chemical_shift_res_list 2-7 10-15

###UUAC_tetraloop
#~/SWA_RNA_python/SWA_dagman_python/SWA_DAG/SWA_rna_minimize.py -silent_file WITH_SHIFT_STATS_Feb_06_SWA.out     -force_field_file rna/rna_12X_CS_03022011_w_intra_terms.wts    -chemical_shift_file renumbered_AUTO_UUAC_shift.bmrb   -chemical_shift_res_list DEFAULT  -common_args common_args_region_FINAL.out -native_pdb G1U4C8_UCAC_1s72.pdb  -num_slave_nodes 500    -LOCATION BIOX2 -remove_virtual_res_variant_during_minimize true > SWA_rna_minimize_LOG.txt 

#Feb_06_WITH_SHIFT_STATS_rebuild_bulge_region_FINAL.out adjusted_AUTO_UUAC_shift.bmrb                          lower_VDW_rep_screener.pdb
#G1U4C8_UCAC_1s72.pdb                                   common_args_region_FINAL.out                           lower_elements.pdb
#WITH_SHIFT_STATS_Feb_06_SWA.out                        fasta                                                  renumbered_AUTO_UUAC_shift.bmrb

#G1U4C8_UCAC_1s72.pdb                                                      adjusted_AUTO_UUAC_shift.bmrb                                             helix_stub_0.pdb
#Nov_10_WITH_SHIFT_STATS_rebuild_bulge_clustered_FINAL_filtered_energy.out common_args_region_FINAL.out                                              helix_stub_0_VDW_rep_screener.pdb
#WITH_SHIFT_STATS_Nov_10_FARFAR.out                                        helix_stub_0.out                                                          renumbered_AUTO_UUAC_shift.bmrb


###LOCAL###:
#SWA_rna_minimize.py -silent_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/WITH_SHIFT_STATS_Jan_23_SWA.out  -tag_range 0 5  -force_field_file rna/rna_M1X_CS_03022011_w_intra_terms.wts     -chemical_shift_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/renumbered_short_bmr5705.bmrb   -chemical_shift_res_list 2-7  -common_args common_args_region_FINAL.out -native_pdb short_2koc_RNA_A.pdb 

#-force_field_file stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts 
#-force_field_file rna/rna_12X_CS_03022011_w_intra_terms.wts 
#-force_field_file rna/rna_4X_CS_03022011_w_intra_terms.wts 
#-force_field_file rna/rna_1X_CS_03022011_w_intra_terms.wts 
#-force_field_file rna/rna_M1X_CS_03022011_w_intra_terms.wts 
#-force_field_file rna/rna_infinitesimal_CS_03022011_w_intra_terms.wts


#-deriv_check false

#########################################################

copy_argv=copy.deepcopy(argv)

#########################################################
chemical_shift_file=parse_options( argv, "chemical_shift_file", "" ) #BMRB_file

if(chemical_shift_file!=""): 
	if(exists(chemical_shift_file)==False): error_exit_with_message("chemical_shift_file (%s) doesn't exist!" %(chemical_shift_file))

#chemical_shift_res_list=parse_seq_num_list_option( argv, "chemical_shift_res_list" )

chemical_shift_res_list=parse_options( argv, "chemical_shift_res_list", [""] ) #DEFAULT or a seq_num_string list.

chemical_shift_H5_prime_mode=parse_options( argv, "chemical_shift_H5_prime_mode", "" )

force_field_file=parse_options( argv, "force_field_file", "" )

in_silent_file=parse_options( argv, "silent_file", "")

in_tag_prefix=parse_options( argv, "in_tag_prefix", "S_")

tag_range=parse_options( argv, "tag_range", [0])

deriv_check=parse_options( argv, "deriv_check", "false") 

perform_o2star_pack=parse_options( argv, "perform_o2star_pack", "true") 

perform_minimize=parse_options( argv, "perform_minimize", "true") #For testing purposes.

rna_torsion_potential_folder=parse_options( argv, "rna_torsion_potential_folder", "")  #ps_03242010/, ps_04282011

remove_virtual_res_variant_during_minimize=parse_options( argv, "remove_virtual_res_variant_during_minimize", "false") 

common_args_file = parse_options(argv, "common_args", "")

native_pdb = parse_options(argv, "native_pdb", "")

LOCATION = parse_options(argv, "LOCATION", "LOCAL")

num_slave_nodes= parse_options(argv, "num_slave_nodes", 500) #FOR LOCATION=BIOX2

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

if(LOCATION not in ["BIOX2", "LOCAL"]): error_exit_with_message("Invalid LOCATION (%s)" %(LOCATION))

if(force_field_file==""): error_exit_with_message("force_field_file==\"\"!")
if(in_silent_file==""): error_exit_with_message("in_silent_file==\"\"!")
if(common_args_file == ""): error_exit_with_message("common_args_file == \"\"!")
if(native_pdb == ""): error_exit_with_message("native_pdb== \"\"!")

if(exists(in_silent_file)==False): error_exit_with_message("in_silent_file (%s) doesn't exist!" %(in_silent_file)) 
if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file)) 
if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!" %(native_pdb) )

##########################################################
#common_args=open( common_args_file  ).readlines()[0][:-1]

common_args=safe_readlines( common_args_file )[0]

common_args=replace_arg_value(common_args, "in:file:silent_struct_type", "binary_rna" ,  allow_missing=True) 

common_args=replace_arg_value(common_args, "output_virtual", "true" ,  allow_missing=True) 

common_args=replace_arg_value(common_args, "silent_read_through_errors", "false" ,  allow_missing=True)

common_args=replace_arg_value(common_args, "constant_seed", "true",  allow_missing=True) 

common_args=replace_arg_value(common_args, "jran", "1111111",  allow_missing=True) 

common_args=replace_arg_value(common_args, "score:weights", force_field_file,  allow_missing=True) 

if(rna_torsion_potential_folder!=""):
	if(use_new_src_code()):
		common_args=replace_arg_value(common_args, "score:rna_torsion_potential", rna_torsion_potential_folder,  allow_missing=True) 
	else:
		common_args=replace_arg_value(common_args, "score:rna_torsion_folder", rna_torsion_potential_folder,  allow_missing=True) 


##########################################################

if(chemical_shift_res_list==""): error_exit_with_message("chemical_shift_res_list==\"\"!")

if(chemical_shift_res_list[0]=="DEFAULT"):

	global_sample_res_list=get_option_name_args_safe(common_args, "global_sample_res_list", [0])

	fasta_file=get_option_name_args_safe(common_args, "fasta", "")

	if(len(global_sample_res_list)==0): error_exit_with_message("global_sample_res_list option does not exist in common_args (%s)" %(common_args))

	if(fasta_file==""): error_exit_with_message("fasta option does not exist in common_rags (%s)" %(common_args))

	sequence = open( fasta_file  ).readlines()[1][:-1]	

	total_res=len(sequence)

	chemical_shift_res_list=add_boundary_seq_num_to_list(global_sample_res_list, total_res, boundary_size=1) 

	
	print_condense_seq_num_list("global_sample_res_list=", global_sample_res_list)
	print "sequence=%s | total_res=%d" %(sequence, total_res)
	print_condense_seq_num_list("Using DEFAULT chemical_shift_res_list=", chemical_shift_res_list)

else:
	chemical_shift_res_list=parse_segment_string_list(chemical_shift_res_list)

	print_condense_seq_num_list("Using USER-SPECIFIED chemical_shift_res_list=", chemical_shift_res_list)

#########################################################
arguments  = ""

arguments += " -algorithm rna_minimize"

arguments += " -database %s" %(get_rosetta_database_folder())

arguments += " -in:file:silent %s " %(in_silent_file)

arguments += " -in:file:native %s " %(native_pdb)

if(perform_o2star_pack=="false"): arguments += " -minimizer_perform_o2star_pack false " 

if(perform_minimize=="false"): 				arguments += " -minimizer_perform_minimize false "

if(deriv_check=="true"):					arguments += " -minimizer_deriv_check true " 

if(remove_virtual_res_variant_during_minimize=="true"): arguments += " -remove_virtual_res_variant_during_minimize true "

if(chemical_shift_file!=""): 
	arguments +=  " -score:rna_chemical_shift_exp_data %s " %(chemical_shift_file)
	if(len(chemical_shift_res_list)>0): 	arguments +=  " -rna_chemical_shift_include_res %s " %(list_to_string(chemical_shift_res_list)) 
	if(chemical_shift_H5_prime_mode!=""): arguments +=  " -score:rna_chemical_shift_H5_prime_mode %s " %(chemical_shift_H5_prime_mode)

arguments += " %s " %(common_args)


if(LOCATION=="BIOX2"):

	if(num_slave_nodes!=0): #If num_slave_nodes=0, means that already created README_SUB.py in the wrapper function!
		create_generic_README_SUB(num_slave_nodes)

	if(exists("CONDOR/")): #This ensures that everything inside this folder is from the latest call to SWA_rna_build_dagman.py
		print "CONDOR/ folder already exist ...removing!"
		submit_subprocess("rm -r CONDOR/")

	system( "mkdir CONDOR/" )  

	if(exists("COMMON_ARGS/")): #This ensures that everything inside this folder is from the latest call to SWA_rna_build_dagman.py
		print "COMMON_ARGS/ folder already exist ...removing!"
		submit_subprocess("rm -r COMMON_ARGS/")

	system( "mkdir COMMON_ARGS/" )  


	fid_dag = open( "rna_build.dag", 'w' ) #Where the dag commands are ouputted.
	fid_dag.write("DOT dag.dot\n")

	minimizer_job_tag="SWA_MINIMIZE"
	minimizer_dag_job_file= '%s/%s.condor' %(DAG_FILE_FOLDER, minimizer_job_tag)

	outfile_dir= "%s/$(Process)/" %(MINIMIZER_FOLDER)
	mapper_outfiles="%s/minimize_silent_file.out " %(outfile_dir)

	arguments += ' -pre_process_output_filename %s ' %( get_rna_minimize_process_output_filename() )
	arguments += ' -out:file:silent %s ' %(mapper_outfiles)
	arguments += ' -job_queue_ID $(Process) ' #determine -tags from job_queue_ID 

	ensure_no_duplicate_options(arguments)

	make_dag_job_submit_file( minimizer_dag_job_file, SWA_RNA_MAIN_EXE, arguments, '', mapper_outfiles , '', '')

	fid_dag.write('\nJOB %s %s\n' % (minimizer_job_tag, minimizer_dag_job_file) )

	swa_minimize_preprocess(in_silent_file, in_tag_prefix, tag_range, minimizer_dag_job_file)

	reducer_outfile="minimize_" + basename(in_silent_file)

	if(len(tag_range)!=0):
		if(len(tag_range)!=2): error_exit_with_message("len(tag_range)!=2")
		reducer_outfile="filtered_range_%d_to_%d_%s" %(tag_range[0], tag_range[1], reducer_outfile) 

	reducer_command =GENERIC_SILENT_FILE_REDUCER_SCRIPT
	reducer_command+=" -reducer_outfile %s -outfolder %s -condor_submit_file %s " %(reducer_outfile, MINIMIZER_FOLDER, minimizer_dag_job_file)

	fid_dag.write('SCRIPT POST %s %s \n' %(minimizer_job_tag, reducer_command) )

	fid_dag.close()

	minimizer_dag_job_lines=open( minimizer_dag_job_file ).readlines()
	
	print 
	print "-------------------%s--------------------" %(minimizer_dag_job_file)
	for line in minimizer_dag_job_lines:
		print line,
	print "-------------------%s--------------------" %(minimizer_dag_job_file)
	print 

elif(LOCATION=="LOCAL"): 

	if(len(tag_range)!=2): error_exit_with_message("len(tag_range)!=2")

	start_tag_num=tag_range[0]
	end_tag_num=tag_range[1]

	if(start_tag_num>end_tag_num): error_exit_with_message("start_tag_num>end_tag_num")
	if(start_tag_num<0): error_exit_with_message("start_tag_num<0")

	arguments += " -VERBOSE true "
	arguments += " -in:file:tags %s %s %s" %(in_tag_prefix, start_tag_num, end_tag_num)

	out_silent_file="filtered_range_%d_to_%d_minimize_%s" (tag_range[0], tag_range[1], basename(in_silent_file))

	arguments += " -out:file:silent %s " %(out_silent_file)

	ensure_no_duplicate_options(arguments)

	command = SWA_RNA_MAIN_EXE + arguments
	command += " > swa_rna_minimize_local.out 2> swa_rna_minimize_local.err"

	print command
	submit_subprocess( command )

else:
	error_exit_with_message("Invalid LOCATION (%s)" %(LOCATION))


print "----------------------------------------------------------------------------------------------------------------------------"
print "----------------------Successfully RAN: %s----------------------" %(list_to_string(copy_argv))
print "----------------------------------------------------------------------------------------------------------------------------"



