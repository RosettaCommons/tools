#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
import math

from SWA_trace_pathway_util import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options

##############Feb 20 2012 Flash Slides#############################################################################

#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_STATS_score_max_n_struct_1000_rebuild_bulge_region_FINAL.out  -rmsd_colname O_loop_rmsd  -y_max 60 -y_min N100 -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_STATS_rebuild_bulge_region_FINAL.out  WITH_SHIFT_STATS_score_max_n_struct_1000_rebuild_bulge_region_FINAL.out  -rmsd_colname O_loop_rmsd  -y_max 60 -y_min N100 -x_max 8.0 

##############Jan 30 2012 Flash Slides#############################################################################

#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_NO_SHIFT_STATS_combine.out  -rmsd_colname NEW_Full_L_rmsd  -y_max N85 -y_min N94  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_Full_L_rmsd  -y_max N40 -y_min N70 -x_max 8.0 


#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -x_label 'RMSD (\\305) to Top-Ranked Rosetta Model' -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 

##############Jan 16 2012 Group Meeting#############################################################################

####UCAC FOX Tetraloop#######

#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_O_loop_rmsd -score_colname shift_RMSD -x_label 'RMSD (\\305) to Top-Ranked Rosetta Model' -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 


####UUAAGU_Hexaloop/######### (A little hacky, use  NEW_NAT_rmsd instead of the usual NEW_O_loop_rmsd)
#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_NO_SHIFT_combine.out  -rmsd_colname NEW_O_loop_rmsd  -y_max N52 -y_min N70  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname  NEW_NAT_rmsd -y_max 80 -y_min N30  -x_max 8.0


#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname  NEW_NAT_rmsd -y_max 30 -y_min N30  -x_max 8.0

#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_NAT_rmsd -score_colname shift_RMSD -x_label 'RMSD (\\305) to Top-Ranked Rosetta Model' -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 

#####TRNA_PHE################

#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_NO_SHIFT_combine.out  -rmsd_colname NEW_O_loop_rmsd  -y_max N60 -y_min N72  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_O_loop_rmsd  -y_max N0 -y_min N30  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_NAT_rmsd -score_colname shift_RMSD -x_label 'RMSD (\\305) to Top-Ranked Rosetta Model' -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 


#####HIV_TAR_LOOP################

#create_score_vs_rmsd_plot.py -silent_file NO_shift_combine.out  -rmsd_colname O_loop_rmsd  -y_max N53.5 -y_min N66  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_STATS_combine.out  -rmsd_colname O_loop_rmsd  -y_max N0 -y_min N40  -x_max 8.0 


#create_score_vs_rmsd_plot.py -silent_file  recal_rmsd_WITH_SHIFT_STATS_combine.out  -rmsd_colname NEW_NAT_rmsd -score_colname shift_RMSD -x_label 'RMSD (\\305) to Top-Ranked Rosetta Model' -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.35 -y_min 0.15  -x_max 8.0 





#####################################################################################################################

###/Dec_19_COMBINE_TRNA_ANTICODON/BLIND_TRNA_A_PLUS_BS_UCC/TRAIL_4_CURR_CODEBASE_CONSISTENCY_CHECK/###
#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.35 -y_min 0.16  -x_max 8.0 


###UGAC (WITHOUT CLOSING WC)
#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 

###WITH_WC_UCAC
#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 

#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.625 -y_min 0.225  -x_max 8.0 


#create_score_vs_rmsd_plot.py -silent_file recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 

##USE THIS##
#create_score_vs_rmsd_plot.py -silent_file suite_clustering_cluster_70_recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.40 -y_min 0.15  -x_max 8.0 



#create_score_vs_rmsd_plot.py -silent_file suite_clustering_cluster_70_recal_rmsd_WITH_SHIFT_STATS_combine.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.25 -y_min 0.15  -x_max 8.0 


#create_score_vs_rmsd_plot.py -silent_file recal_full_WITH_SHIFT_STATS_region_FINAL.out -rmsd_colname NEW_Full_L_rmsd -score_colname shift_RMSD -y_label 'Chemical Shift RMSD (ppm)' -y_max 0.6 -y_min 0.2  -x_max 8.0 


#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_RMSD_cluster_75_FINAL_filtered_energy.out -y_max N30 -y_min N60  -x_max 8.0
#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_RMSD_region_FINAL.out  -y_max N30 -y_min N60  -x_max 8.0

#create_score_vs_rmsd_plot.py -silent_file cluster_75_FINAL_filtered_energy.out -y_max N55 -y_min N67  -x_max 8.0
#create_score_vs_rmsd_plot.py -silent_file region_FINAL.out -y_max N55 -y_min N67  -x_max 8.0

#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_RMSD_score_max_n_struct_2500_cluster_0_FINAL_filtered_energy.out -y_max N140.0 -y_min N154.0  -x_max 6.0

#create_score_vs_rmsd_plot.py -silent_file score_max_n_struct_2500_cluster_0_FINAL_filtered_energy.out -y_max N150.0 -y_min N160.0  -x_max 6.0 

#create_score_vs_rmsd_plot.py -silent_file region_FINAL.out   -y_max N106.0 -y_min N120.0 -x_max 6.0 

#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_RMSD_region_FINAL.out   -y_max N80.0 -y_min N112.0 -x_max 6.0 

#create_score_vs_rmsd_plot.py -silent_file region_FINAL.out suite_clustering_cluster_70_FINAL_filtered_energy.out  -y_max N40.0 -y_min N49.0 -x_max 6.0 

#create_score_vs_rmsd_plot.py -silent_file WITH_SHIFT_RMSD_region_FINAL.out  WITH_SHIFT_RMSD_suite_clustering_cluster_70_FINAL_filtered_energy.out -y_max N20.0 -y_min N37.0 -x_max 6.0 

#####################################################



rmsd_colname=parse_options(argv, "rmsd_colname", "O_loop_rmsd")
score_colname=parse_options(argv, "score_colname", "score")

silent_file_list=parse_options(argv, "silent_file", [""])
point_color=parse_options(argv, "point_color", "red")
point_size=parse_options(argv, "point_size", 1.5)

x_label=parse_options(argv, "x_label", "RMSD (\\305)")
y_label=parse_options(argv, "y_label", "Energy Score")

#y_label=list_to_string(y_label)[1:]

user_y_min=parse_options(argv, "y_min", "NONE")
user_y_max=parse_options(argv, "y_max", "NONE")
user_x_min=parse_options(argv, "x_min", "0.0")
user_x_max=parse_options(argv, "x_max", "8.0")

if(len(argv)!=1): error_exit_with_message("len(argv)!=1, leftover_argv=%s" %(list_to_string(argv) ) )

for silent_file in silent_file_list:
	if(exists(silent_file)==False): error_exit_with_message("silent_file (%s) doesn't exist!" %(silent_file))

	assert_is_valid_non_empty_silent_file( silent_file )


figure_name="score_vs_rmsd_%s.eps" %(basename(silent_file)) 

print "figure_name= %s" %(figure_name)

script_name='gnuplot_script.txt'

system( 'rm %s' % script_name)  #create parent directory as needed
system( 'rm %s' % figure_name)  #create parent directory as needed

outfile = open( script_name, 'w' ) 


outfile.write('\nset terminal postscript enhanced color\n')
outfile.write('set output "%s"\n' % figure_name)


point_type_ID=7

if(point_color=="red"):
	outfile.write('set style line 7 lt 1 lw 2\n') #line 7-> filled circle lt 1 red
elif(point_color=="blue"):
	outfile.write('set style line 7 lt 3 lw 2\n') #line 7-> filled circle lt 3 blue
else:
	error_exit_with_message("Invalid point_color (%s)" %(point_color))

'''
if(figure_ID==2 or figure_ID==13): 
	outfile.write('set style line 6 lt 3 lw 2\n')	#line 6-> hollow circle lt 3 blue 
	outfile.write('set style line 7 lt 1 lw 2\n') #line 7-> filled circle lt 1 red
	point_type_ID=6
elif(figure_ID==9):
	outfile.write('set style line 7 lt 3 lw 2\n') #line 7-> filled circle lt 3 blue
	outfile.write('set style line 6 lt 1 lw 6\n') #line 6-> hollow circle lt 1 red
elif(figure_ID==11):
	outfile.write('set style line 7 lt 3 lw 2\n') #line 7-> filled circle lt 3 blue
	outfile.write('set style line 31 lt 1 lw 0.0\n') #line 31-> filled circle lt 1 red
	point_type_ID=7
elif(figure_ID==12):
	outfile.write('set style line 7 lt 4 lw 2.0\n') #line 6-> hollow circle lt 1 red
	outfile.write('set style line 31 lt 1 lw 0.0\n') #line 31-> filled circle lt 1 red
	point_type_ID=7
else: #standard
	outfile.write('set style line 7 lt 3 lw 2\n') #line 7-> filled circle lt 3 blue
	outfile.write('set style line 6 lt 1 lw 2\n') #line 6-> hollow circle lt 1 red
'''


outfile.write('set nokey\n')
outfile.write('set bmargin at screen %8.4f\n' %(0.14) )
outfile.write('set size %f,%f\n' %(1.0, 1.0) )  
outfile.write('set pointsize %f\n' %(point_size))


outfile.write('set ylabel font "Helvetica,30"\n')
outfile.write('set xlabel font "Helvetica,30"\n')

outfile.write('set encoding iso_8859_1\n')

outfile.write('set ylabel "%s"\n' %(y_label)) #April 7, 2011
outfile.write('set xlabel "%s"\n' %(x_label))

outfile.write('set ylabel offset -4.5,0.0\n')
outfile.write('set xlabel offset 0.0,-1.0\n')

outfile.write('set xtics font "Helvetica,30"\n')
outfile.write('set ytics font "Helvetica,30"\n')
outfile.write('set xtics 1.0\n')


global_y_min=10000000000000.000
global_y_max=-10000000000000.000

plot_command='plot ' 

for n in range(len(silent_file)):

	silent_file=silent_file_list[n]

	rmsd_column=get_silent_file_col_index(rmsd_colname, silent_file)+1
	score_column=get_silent_file_col_index(score_colname, silent_file)+1

	(y_min,y_max)=get_min_max(silent_file) #TIME CONSUMING...
	
	if(global_y_min>y_min): global_y_min=math.floor(y_min)

	if(global_y_max<y_max): global_y_max=math.floor(y_max)+1

	plot_command+=' \'<grep "SCORE" %s\'	using %d:%d ls %d,' % (silent_file ,rmsd_column, score_column, point_type_ID) 

plot_command=plot_command[:-1]


if(user_y_min!="NONE"):
	global_y_min=convert_string_to_float(user_y_min)

if(user_y_max!="NONE"):
	global_y_max=convert_string_to_float(user_y_max)

x_min=convert_string_to_float(user_x_min)
x_max=convert_string_to_float(user_x_max)


outfile.write('set yrange[%f:%f]\n' %(global_y_min, global_y_max))
outfile.write('set xrange[%f:%f]\n' %(x_min, x_max))

outfile.write(plot_command + "\n")

outfile.write('set output\n')
outfile.write('reset\n\n')

outfile.close()

print "PLOTTING.....%s" %(figure_name) 
sys.stdout.flush()

submit_subprocess('gnuplot %s' %(script_name))
submit_subprocess('gv %s &' %(figure_name) )

print "---------------------START OF %s----------------------------" %(script_name)
submit_subprocess('more %s & ' %(script_name))
sleep(0.05)
print "---------------------END OF   %s----------------------------" %(script_name)


