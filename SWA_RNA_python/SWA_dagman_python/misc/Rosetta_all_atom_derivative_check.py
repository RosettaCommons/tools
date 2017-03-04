#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options, parse_seq_num_list_option
######################################################################


#plot '<grep "syn_chi_torsion"  output.txt ' u 2:4
#syn_chi_torsion= -25.4147 syn_chi_offset= 24.5853

#plot '<grep "ratio   5"   rosetta_deriv_check_log.out  ' u 11:12, x
#plot '<grep "ratio   5"   rosetta_deriv_check_log.out  ' u 11:($12-$11), 0
#plot '<grep "ratio   5"   rosetta_deriv_check_log.out  ' u 11:(100*($12-$11)/$12), 100, -100, 0 ls 3

#core.optimization: ratio inc rsd typ atm prsd ptyp patm natm   numeric  analytic     ratio       f11       f00       f22  vars[ii]
#numeric--> column 11
#analytic--> column 12

#########################################################

#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/Rosetta_rna_file/SWA_LOOP_PAPER/RELEASE_EXAMPLE_1S72_23_rRNA_35_40/native.pdb -force_field_file rna/test.wts  -rna_torsion_potential_folder ps_03242010/


#CD_geom_sol_intra_RNA_rna_hires_07232011_with_intra_base_phosphate.wts
#rna_hires_07232011_with_intra_base_phosphate.wts 

###BIOX2:
#bsub -n 1  ~/SWA_RNA_python/SWA_dagman_python/misc/Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/BIOX2_Feb_06_compare_CD_and_CI_geom_sol_intra_RNA/short_TL1_1f7y_RNA.pdb  -force_field_file rna/test.wts  -min_type dfpmin  -deriv_check False 

###ADE:
#~/SWA_RNA_python/SWA_dagman_python/misc/Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/ADE_Feb_06_compare_CD_and_CI_geom_sol_intra_RNA/short_TL1_1f7y_RNA.pdb  -force_field_file rna/test.wts  -min_type dfpmin  -deriv_check False


###LOCAL:
#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/short_2koc_RNA_A.pdb   -force_field_file rna/rna_CS_03022011_w_intra_terms.wts    -chemical_shift_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/renumbered_short_bmr5705.bmrb  -deriv_check False -chemical_shift_res_list 2-7  -min_type dfpmin

#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/short_2koc_RNA_A.pdb   -force_field_file stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts     -chemical_shift_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/renumbered_short_bmr5705.bmrb  -deriv_check False -chemical_shift_res_list 2-7  -min_type dfpmin


#Rosetta_all_atom_derivative_check.py -silent_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/WITH_SHIFT_STATS_Jan_23_SWA.out  -silent_tag S_1413 -force_field_file rna/rna_12X_CS_03022011_w_intra_terms.wts          -chemical_shift_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/renumbered_short_bmr5705.bmrb  -deriv_check False -chemical_shift_res_list 2-7  -min_type dfpmin

#-force_field_file rna/rna_CS_03022011_w_intra_terms.wts  
#-force_field_file stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts  



#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/short_2koc_RNA_A.pdb   -force_field_file rna/chem_shift_only.wts   -chemical_shift_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/renumbered_short_bmr5705.bmrb  -deriv_check True -chemical_shift_res_list 2-7  

#-chemical_shift_H5_prime_mode UNIQUE -min_type dfpmin

#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/short_2koc_RNA_A.pdb   -force_field_file rna/fa_stack_only.wts  -chemical_shift_file ~/minirosetta/test/Feb_22_Rosetta_chem_shift_CODE/renumbered_short_bmr5705.bmrb -chemical_shift_H5_prime_mode UNIQUE -deriv_check True -chemical_shift_res_list 2-7 







#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/Feb_06_compare_CD_and_CI_geom_sol_intra_RNA/short_TL1_1f7y_RNA.pdb  -force_field_file rna/test.wts  -min_type dfpmin  -deriv_check False

#Rosetta_all_atom_derivative_check.py -silent_file ~/minirosetta/test/Feb_04_COMPARE_TORSION_WITHOUT_AND_WITH_rna_torsion_skip_chainbreak/Jan_07_2012_SWA_ss_loop.out  -force_field_file stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts -rna_torsion_potential_folder ps_04282011  -deriv_check False -rna_torsion_skip_chainbreak true

#rna_torsion_skip_chainbreak			Oct_20_2011_SWA.out   Jan_28_FARFAR.out 

#Rosetta_all_atom_derivative_check.py -pdb ~/minirosetta/test/Feb_04_COMPARE_OLD_AND_NEW_SWA_IDEALIZE_HELIX_ENERGY/RESETTED_BEFORE_o2star_pack_OLD_SWA_helix.pdb  -force_field_file stepwise/rna/rna_hires_07232011_with_intra_base_phosphate.wts   -rna_torsion_potential_folder ps_03242010/ -deriv_check False 

#OLD_SWA_helix.pdb   NEW_SWA_helix.pdb RESETTED_BEFORE_o2star_pack_NEW_SWA_helix.pdb  RESETTED_BEFORE_o2star_pack_OLD_SWA_helix.pdb

#Rosetta_all_atom_derivative_check.py -pdb /Users/sripakpa/minirosetta/test/01_2012_Deriv_test/short_TL1_1f7y_RNA.pdb -force_field_file rna/test_2.wts   -rna_torsion_potential_folder ps_03242010/ -deriv_check False 

#Rosetta_all_atom_derivative_check.py -pdb /Users/sripakpa/minirosetta/test/01_2012_Deriv_test/short_TL1_1f7y_RNA.pdb -force_field_file rna/test_2.wts   -rna_torsion_potential_folder ps_03242010/ -deriv_check False 

#Rosetta_all_atom_derivative_check.py -silent_file /Users/sripakpa/minirosetta/test/01_2012_Deriv_test/CURR_CODE_SWA_1S72.out -force_field_file rna/test_2.wts   -rna_torsion_potential_folder ps_03242010/ -perform_minimizer_run false -skip_o2star_trials false

#Rosetta_all_atom_derivative_check.py -silent_file /Users/sripakpa/minirosetta/test/01_2012_Deriv_test/CURR_CODE_SWA_1S72.out -force_field_file rna/test_2.wts   -rna_torsion_potential_folder ps_03242010/ -deriv_check False 

#Rosetta_all_atom_derivative_check.py -pdb /Users/sripakpa/minirosetta/test/01_2012_Deriv_test/1S72_loop_native.pdb -force_field_file rna/test.wts  -rna_torsion_potential_folder ps_03242010/

#Rosetta_all_atom_derivative_check.py -pdb /Users/sripakpa/minirosetta/test/01_2012_Deriv_test/short_TL1_1f7y_RNA.pdb -force_field_file rna/test.wts -min_type dfpmin  -deriv_check False 



#########################################################

copy_argv=copy.deepcopy(argv)

#########################################################
chemical_shift_file=parse_options( argv, "chemical_shift_file", "" ) #BMRB_file

chemical_shift_res_list=parse_seq_num_list_option( argv, "chemical_shift_res_list" )

chemical_shift_H5_prime_mode=parse_options( argv, "chemical_shift_H5_prime_mode", "" )

force_field_file=parse_options( argv, "force_field_file", "" )

pdb_file=parse_options( argv, "pdb", "") 

silent_file=parse_options( argv, "silent_file", "")

silent_tag=parse_options( argv, "silent_tag", "")

deriv_check=parse_options( argv, "deriv_check", "True") 

min_type=parse_options( argv, "min_type", "linmin") #"dfpmin"

skip_o2star_trials=parse_options( argv, "skip_o2star_trials", "true") #Avoid randomness

perform_minimizer_run=parse_options( argv, "perform_minimizer_run", "true") #For testing purposes.

rna_torsion_skip_chainbreak=parse_options( argv, "rna_torsion_skip_chainbreak", "true" ) #Implemented code on Jan 19, 2012; Set to false to check backward compatibility.

rna_torsion_potential_folder=parse_options( argv, "rna_torsion_potential_folder", "")  #ps_03242010/, ps_04282011

#-score:rna_torsion_potential ps_03242010/

if( force_field_file==""): error_exit_with_message('force_field_file==""')

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )
##########################################################

#if(OLD_CODE): exe_folder=HOMEDIR + "/src/March_24_2011_BACKUP_COPY_mini/bin/"


command=""
if(use_new_src_code()):
	command+= get_rosetta_EXE("swa_rna_util") 

else:
	command+= get_rosetta_EXE("parin_test") 
	

database_folder= get_rosetta_database_folder()

#########################################################
command+= " -algorithm rna_fullatom_minimize_test"
command += " -database %s" %(database_folder)

num_struct_src=0

if(pdb_file!=""):
	num_struct_src+=1
	if( exists( pdb_file )==False):  error_exit_with_message("pdb_file (%s) doesn't exist!" %(pdb_file) )
	command += " -input_tag_list " + pdb_file

if(silent_file!=""):
	num_struct_src+=1
	if( exists( silent_file )==False):  error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file) )
	command += " -in:file:silent %s " %(silent_file)

	if(silent_tag!=""): command += " -tags %s " %(silent_tag)
	

if(num_struct_src!=1): error_exit_with_message("num_struct_src!=1")

command += " -score:weights %s " %(force_field_file) 
command += " -minimizer_min_type %s " %(min_type)

if(rna_torsion_potential_folder!=""):

	if(use_new_src_code()):
		command += " -score:rna_torsion_potential %s " %(rna_torsion_potential_folder)
	else:
		command += " -score:rna_torsion_folder %s " %(rna_torsion_potential_folder)

command += " -minimizer_skip_o2star_trials %s " %(skip_o2star_trials)

command += " -minimizer_perform_minimizer_run %s " %(perform_minimizer_run)

if(deriv_check==False): command+= " -minimizer_deriv_check false "

if(chemical_shift_file!=""): 
	if(exists(chemical_shift_file)==False): error_exit_with_message("chemical_shift_file (%s) doesn't exist!" %(chemical_shift_file))
	command+= " -score:rna_chemical_shift_exp_data %s " %(chemical_shift_file)
	if(len(chemical_shift_res_list)>0): command+= " -rna_chemical_shift_include_res %s " %(list_to_string(chemical_shift_res_list)) 

if(chemical_shift_H5_prime_mode!=""): command+= " -score:rna_chemical_shift_H5_prime_mode %s " %(chemical_shift_H5_prime_mode)

command+= " -score:rna_torsion_skip_chainbreak %s " %(rna_torsion_skip_chainbreak) #Feb 04, 2012.

command += " -output_virtual true " #Added this on Jan 19, 2012

command += " -constant_seed true -jran 1111111 " #Added this on Jan 20, 2012..This 'should' prevent randomness in the o2star_packer!

command += " > rosetta_deriv_check_log.out 2> rosetta_deriv_check_log.err"
print command
submit_subprocess( command )

print "----------------------------------------------------------------------------------------------------------------------------"
print "----------------------%s----------------------" %(list_to_string(copy_argv))
print "Rosetta_all_atom_derivative_check.py successfully RAN! "
print "----------------------------------------------------------------------------------------------------------------------------"








































#########################################################

'''
						ratio_header_output = true;
2						TR << "ratio" << 
3							A( 4, "inc" ) <<
4							A( 4, "rsd" ) <<
5							A( 4, "typ" ) <<
6							A( 4, "atm" ) <<
7							A( 5, "prsd" ) <<
8							A( 5, "ptyp" ) <<
9							A( 5, "patm" ) <<
10							A( 5, "natm" ) <<
11							A( 10, "numeric" ) <<
12							A( 10, "analytic" ) <<
13							A( 10, "ratio" ) <<
14							A( 10, "f11" ) <<
15							A( 10, "f00" ) <<
16							A( 10, "f22" ) <<
17							A( 10, "vars[ii]" ) << std::endl;
					}


2					TR << "ratio" <<
3						I( 4, j ) <<
4						I( 4, dof_node.rsd() ) <<
5						I( 4, dof_node.type() ) <<
6						I( 4, dof_node.atomno() ) <<
7						I( 5, parent_id.rsd() ) <<
8						I( 5, parent_id.type() ) <<
9						I( 5, parent_id.atomno() ) <<
10						I( 5, dof_node.atoms().size()) <<
11						F( 10, 4, deriv ) <<                // column 11
12						F( 10, 4, dE_dvars[ii] ) <<         // column 12
13						F( 10, 4, ratio ) <<
14						F( 10, 4, f11 ) <<
15						F( 10, 4, f00 ) <<
16						F( 10, 4, f22 ) <<
17						F( 10, 4, start_vars[ii] ) << std::endl;
				}
'''



#########################################################
'''
sripakpa@rescomp-09-171908:~/minirosetta$ find . -name 'syn_chi_pose.pdb'
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_13_DERIV_CHECK_GEOM_SOLV/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_14_DERIV_CHECK_FA_INTRA_RNA/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_16_DERIV_CHECK_GEOM_SOLV_VERBOSE/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_17_DERIV_CHECK_GEOM_SOLV_VERBOSE/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_18_DERIV_CHECK_GEOM_SOLV_VERBOSE_EXCLUDE_O3STAR/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_19_DERIV_CHECK_GEOM_SOLV_VERBOSE_EXCLUDE_O3STAR/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_20/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_21/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_23/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_26_HBOND_DERIV_CHECK/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_27/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_28_SMALLER_INCREMENT_0.000005/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_29/syn_chi_pose.pdb
./Rosetta_rna_file/TURNER_GAGU_DUPLEX/06_2011/June_29_DEBUG_intra_base_phospate/test_1/TRAIL_31_SMALLER_INCREMENT_0.000005_EVERY_TERMS/syn_chi_pose.pdb
'''
####OLD####
#Rosetta_all_atom_derivative_check.py -pdb expand_radius_50_3d2v_RNA_A.pdb -force_field_file test.wts

#Rosetta_all_atom_derivative_check.py -pdb short_TL1_1f7y_RNA.pdb -force_field_file torsion_only.wts

#Rosetta_all_atom_derivative_check.py -pdb syn_chi_pose.pdb -force_field_file rna_hires_06262011_with_intra_base_phosphate.wts
####OLD####
