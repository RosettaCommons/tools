#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os
import copy
######################################################################

from SWA_util import *
from SWA_parse_options import parse_options
######################################################################


# plot '<grep "new_cluster"  cluster_rotamers_RMSD_3.0.out ' using 3:5


#cluster_rotamers.py  -optimize_screening True -cluster_rmsd 1.0 -bin_size 10 -two_stage_clustering True

#cluster_rotamers.py  -optimize_screening True -cluster_rmsd 1.5 -two_stage_clustering False   -VDW_rep_screen_slow_check true 

# plot '<grep "new_cluster"  TRAIL_3_OPTIMIZE/cluster_rotamers_RMSD_2.0.out ' using 7:3, '<grep "new_cluster"  TRAIL_4_OPTIMIZE_1_ANGSTROM/cluster_rotamers_RMSD_1.0.out ' using 7:3

copy_argv=copy.deepcopy(argv)


rotamer_cluster_rmsd=parse_options( argv, "cluster_rmsd", "0.0" ) 
graphic=parse_options( argv, "graphic", "false" ) 
create_silent_file=parse_options( argv, "create_silent_file", "False" ) 
analyze_silent_file=parse_options( argv, "analyze_silent_file", "False" ) 
sequence=parse_options( argv, "sequence", "aa" ) 
cluster_rotamers_optimize_screening=parse_options( argv, "optimize_screening", "True")
two_stage_clustering=parse_options( argv, "two_stage_clustering", "True")
quick_test=parse_options( argv, "quick_test", "false")
bin_size=parse_options( argv, "bin_size", 20)
VDW_rep_screen=parse_options( argv, "VDW_rep_screen", "true")
VDW_rep_screen_slow_check=parse_options( argv, "VDW_rep_screen_slow_check", "false")

if(rotamer_cluster_rmsd<0.01): error_exit_with_message("rotamer_cluster_rmsd=(%s)<0.01" %(rotamer_cluster_rmsd) )

if(create_silent_file): two_stage_clustering=False

if(len(argv)!=1):
	print argv, " leftover len(argv)=", len(argv)
	assert(False)


#####################################

EXE = get_rosetta_EXE("parin_test") 

database_folder= get_rosetta_database_folder()

####################################

output_silent_file="rotamer_silent_file_%s.out" %(sequence)

if(analyze_silent_file):
	print "analyze the silent_file"

else:
	command=EXE
	command += " -algorithm cluster_rotamers"
	command += " -database %s" %(database_folder)
	command += " -rotamer_cluster_rmsd %s " %(rotamer_cluster_rmsd)
	command += " -output_virtual false "
	command += " -graphic %s " %(graphic)
	command += " -dinucleotide_sequence %s " %(sequence)
	command += " -cluster_rotamers_optimize_screening %s " %(cluster_rotamers_optimize_screening)
	command += " -two_stage_rotamer_clustering %s " %(two_stage_clustering)
	command += " -quick_test %s " %(quick_test)
	command += " -cluster_rotamer_bin_size %s " %(bin_size)
	command += " -cluster_rotamer_replusion_screen %s " %(VDW_rep_screen)
	command += " -cluster_rotamer_VDW_rep_screening_slow_check %s " %(VDW_rep_screen_slow_check)

	if(float(rotamer_cluster_rmsd)<1.01): command += " -cluster_rotamer_sparse_output true "

	if(create_silent_file): 
		if(exists(output_silent_file)): 
			print "output_silent_file (%s) already exist! ...removing.." %(output_silent_file)
			submit_subprocess("rm %s " %(output_silent_file) )
		command+= " -output_silent_file %s " %(output_silent_file)

	command += " > cluster_rotamers_RMSD_%s.out 2> cluster_rotamers_RMSD_%s.err" %(rotamer_cluster_rmsd, rotamer_cluster_rmsd)
	print command
	submit_subprocess( command )



print "----------------------------------------------------------------------------------------------------------------------------"
print "----------------------%s----------------------" %(list_to_string(copy_argv))
print "cluster_rotamers.py sucessfully RAN! "
print "----------------------------------------------------------------------------------------------------------------------------"


'''
Used c++ code:

		/*Problem is that clash between base and sugar is not well defined...I think the only possible clash is between the phosphate and the (sugar+base).
		//Should just manually write this up..
		//maybe should set anchor nuceltoide to A-form and screen for clashes.
		if(replusion_screen){
			(*scorefxn)(pose);

			EnergyMap const & energy_map = pose.energies().total_energies();

			Real const rep_score = scorefxn->get_weight(fa_rep) * energy_map[scoring::fa_rep];
			Real const intra_rep_score = scorefxn->get_weight(fa_intra_rep) * energy_map[scoring::fa_intra_rep];
			
			if((intra_rep_score+rep_score)>(1.0/0.12)) continue;
			pass_rep_screen_count++;
			
		}
		*/
'''
