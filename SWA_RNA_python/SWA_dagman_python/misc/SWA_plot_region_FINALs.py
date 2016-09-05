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
######################################################################

from SWA_trace_pathway_util import *
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################

#A-A platofrm
#SWA_plot_region_FINALs.py -quick False  -recal_rmsd False -num_elements 3  > Sept_6_SWA_plot_region_FINAL_TRAIL.log  2> Sept_6_SWA_plot_region_FINAL_TRAIL.err



#J5/Ja
#SWA_plot_region_FINALs.py -quick False -recal_rmsd True -num_elements 6  > April_8_SWA_plot_region_FINAL_TRAIL.log  2> April_8_SWA_plot_region_FINAL_TRAIL.err

#3D2V
#SWA_plot_region_FINALs.py -quick False -recal_rmsd True -num_elements 6  > plot_region_FINAL_TRAIL.log  2> plot_region_FINAL_TRAIL.err


#3D2V //Farfar + SWA_native_cst
#SWA_plot_region_FINALs.py -quick False -recal_rmsd False -num_elements 6  > plot_region_FINAL_TRAIL.log  2> plot_region_FINAL_TRAIL.err


#SWA_plot_region_FINALs.py -quick False -recal_rmsd False -num_elements 6 

num_elements= parse_options( argv, "num_elements", 0)
test_mode= parse_options( argv, "test_mode", "False")
one_point_type= parse_options( argv, "one_point_type", "True")
quick= parse_options( argv, "quick", "True")
RAW_SCORE_FILE= parse_options( argv, "RAW_SCORE_FILE", "True")
recal_rmsd= parse_options( argv, "recal_rmsd", "False")
user_input_rmsd_colname= parse_options( argv, "rmsd_colname", "")
four_edge= parse_options( argv, "four_edge", "False")

###WARNING quick have different meaning depending on whether the recal_rmsd silent_files or the original silent_files are used
###In fact (recal_rmsd=true + quick=false) should contain contain data points as (recal_rmsd=false + quick=true)

rmsd_colname="O_loop_rmsd"
if(recal_rmsd): rmsd_colname="NEW_Full_L_rmsd" ##Used to be "NEW_O_loop_rmsd" before March 28, 2011

if(user_input_rmsd_colname!=""): 
	rmsd_colname=user_input_rmsd_colname
	print "WARNING: user_input_rmsd_colname overwrite DEFAULT rmsd_colname!"


figure_ID=14

if(figure_ID==5 or figure_ID==6): four_edge=True

print "--------------------------------------------------------------------"
print "--------------------------------------------------------------------"
print "--------------------------------------------------------------------"
print "quick= %s" %(quick)
print "recal_rmsd= %s" %(recal_rmsd)
print "rmsd_colname= %s" %(rmsd_colname)
print "four_edge= %s" %(four_edge)
print "--------------------------------------------------------------------"
print "--------------------------------------------------------------------"
print "--------------------------------------------------------------------"




FARFAR_command_before=""
FARFAR_command_after=""

#rmsd_column=27
score_column=2

if(num_elements==0): error_exit_with_message('num_elements==0')  

#SWA_trace_pathway.py -tag  S_0 -silent_file region_FINAL.out

figure_name="FINAL_regions.out_" + os.path.abspath(os.curdir) 
figure_name_inside=figure_name
figure_name_inside=figure_name_inside.replace('_','\_')

figure_name=figure_name.replace('/','_') + '.eps'
print "figure_name= %s" %figure_name


script_name='gnuplot_script.txt'

system( 'rm %s' % script_name)  #create parent directory as needed
system( 'rm %s' % figure_name)  #create parent directory as needed

outfile = open( script_name, 'w' ) 



outfile.write('\nset terminal postscript enhanced color\n')
outfile.write('set output "%s"\n' % figure_name)

point_type_ID=7

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


#outfile.write('set style line 1 linetype "solid" lw  5\n')
#outfile.write('set border lt rgb "blue" lw 1\n')

if(one_point_type): outfile.write('set nokey\n')
#outfile.write('set size 1.0, 1.0\n')
#outfile.write('set origin %f,%f\n' %(0.0, -0.1) )
outfile.write('set bmargin at screen %8.4f\n' %(0.14) )
outfile.write('set size %f,%f\n' %(1.0, 1.0) )  
if(figure_ID==4 or figure_ID==5 or figure_ID==11):
	outfile.write('set pointsize 1.0\n')
else:
	outfile.write('set pointsize 1.5\n')


outfile.write('set ylabel font "Helvetica,30"\n')
outfile.write('set xlabel font "Helvetica,30"\n')

outfile.write('set encoding iso_8859_1\n')

#outfile.write('set ylabel "Energy Score"\n')
outfile.write('set ylabel "Rosetta Energy"\n') #April 7, 2011
outfile.write('set xlabel "RMSD (\\305)"\n')

outfile.write('set ylabel offset -4.5,0.0\n')
outfile.write('set xlabel offset 0.0,-1.0\n')

outfile.write('set xtics font "Helvetica,30"\n')
outfile.write('set ytics font "Helvetica,30"\n')
outfile.write('set xtics 1.0\n')


global_min=10000000000000000
global_max=-10000000000000000


#set multiplot title "Demo of placing multiple plots (2D and 3D)\nwith explicit alignment of plot borders"
plot_command='plot '

file_ID=0

file_name_list=[]

if(four_edge):
	folder_name_list = glob( "REGION_L%d_*" %(num_elements) )	
	print "FOUR EDGE folder_name_list= ", folder_name_list

	for folder_name in folder_name_list:

		if(recal_rmsd):
			error_exit_with_message("recal_rmsd not implemented for four_edge mode yet." )
		else:
			file_name_string="%s/SCORE_start_from_region_l%s_*_sample.out"  %(folder_name, num_elements-1)

			if(quick): 
				file_name_string="%s/start_from_region_l%s_*_sample_filtered.out"  %(folder_name, num_elements-1)
			file_name_list.extend( glob( file_name_string ) )

	print "FOUR EDGE file_name_list= ", file_name_list	


else:
	################################################################################################
	if(quick): 
		filename="region_FINAL.out"
		if(exists(filename)==False): error("filename (%s) doesn't exist!" %(filename) )

		file_name_list=[filename]

		if(RAW_SCORE_FILE): error_exit_with_message("quick==TRUE and RAW_SCORE_FIL==TRUE")

	else:
		for j in range(num_elements):

			#if(len(file_name_list)>0): break

			i=(j+1) % num_elements

			folder_name="REGION_%d_%d/" % (i,j)

			if(exists(folder_name)==False): continue

			print "folder_name= %s" %(folder_name)

			'''
			region_name_list=[]

			for step_i in range(0,3,1):
				for step_j in range(0,3,1):
					filename=""
					region_name="%d_%d" %( (i+step_i) % num_elements , (j-step_j) % num_elements )
			'''


			if(recal_rmsd):
				globstring=folder_name + "/RECAL_RMSD/" + "recal_full_start_from_region_*_sample_filtered.out"  #After March 28, 2011
			else:
				if(RAW_SCORE_FILE):
					globstring=folder_name + "/" + "SCORE_start_from_region_*_sample.out" 
				else:
					globstring=folder_name + "/" + "start_from_region_*_sample_filtered.out" 


			spec_region_file_name_list = glob( globstring )
			spec_region_file_name_list.sort()


			for file_name in spec_region_file_name_list:

				file_name=os.path.abspath(file_name)

				if(exists(file_name)==False): continue 

				#############################################################################	
				if(figure_ID==11):

					SCORE_file_name=dirname(file_name) + "/" + "SCORE_" + basename(file_name)

					if(exists(SCORE_file_name)==False):
						grep_SCORE_command=	'grep "SCORE: " %s > %s' %(file_name, SCORE_file_name)
						#print "grep_SCORE_command= %s" %(grep_SCORE_command)
						sys.stdout.flush()
						submit_subprocess(grep_SCORE_command)

					sleep( 1 )

					energy_cut_file_name=dirname(file_name) + "/" + "SCORE_energy_cut_" + basename(file_name)

					ENERGY_CUT_FILE=open( energy_cut_file_name, 'w' ) 
					
					try:
						infile = open( SCORE_file_name, 'r' )
					except:
						error_exit_with_message("cannot open %s " % SCORE_file_name)

					first_line=True
					line_count=0

					for line in infile:

						#if(line_count>1000): break

						line_count+=1	


						if(Is_column_name_line(line)!=first_line): 
							error_exit_with_message("Is_column_name_line(line)=(%s)!=(%s)=first_line, line=%s" %(Is_column_name_line(line), first_line, line) )
							
						if(first_line):
							first_line=False
							ENERGY_CUT_FILE.write( line )
							continue

						energy_score=float(line.split()[1])

						#SCORE:    -95.107    S_000000

						if(energy_score>(-95.107+14)): continue

						ENERGY_CUT_FILE.write( line )
		
					ENERGY_CUT_FILE.close()
					file_name=energy_cut_file_name

				#############################################################################	
				if(figure_ID==12):

					SCORE_file_name=dirname(file_name) + "/" + "SCORE_" + basename(file_name)

					if(exists(SCORE_file_name)==False):
						grep_SCORE_command=	'grep "SCORE: " %s > %s' %(file_name, SCORE_file_name)
						#print "grep_SCORE_command= %s" %(grep_SCORE_command)
						sys.stdout.flush()
						submit_subprocess(grep_SCORE_command)

					sleep( 1 )

					RMSD_cut_file_name=dirname(file_name) + "/" + "SCORE_RMSD_cut_" + basename(file_name)

					RMSD_CUT_FILE=open( RMSD_cut_file_name, 'w' ) 
					
					try:
						infile = open( SCORE_file_name, 'r' )
					except:
						error_exit_with_message("cannot open %s " % SCORE_file_name)

					first_line=True
					line_count=0

					for line in infile:

						#if(line_count>1000): break

						line_count+=1	


						if(Is_column_name_line(line)!=first_line): 
							error_exit_with_message("Is_column_name_line(line)=(%s)!=(%s)=first_line, line=%s" %(Is_column_name_line(line), first_line, line) )
							
						if(first_line):
							first_line=False
							RMSD_CUT_FILE.write( line )
							continue

						RMSD_value=float(line.split()[26])
						ENERGY_value=float(line.split()[1])

						#SCORE:    -95.107    S_000000

						if(RMSD_value>(1.5)): continue
						if(ENERGY_value>(-91)): continue

						RMSD_CUT_FILE.write( line )
		
					RMSD_CUT_FILE.close()
					file_name=RMSD_cut_file_name
				#############################################################################

				print "file_name= %s" %(file_name)
				sys.stdout.flush()
				file_name_list.append(file_name)
				

	################################################################################################


################################################################################################
file_name_count=0
if(test_mode): file_name_list=[file_name_list[0]]

print "------------------------------------------------------------------------------"
print "file_name_list="
for file_name in file_name_list:
	file_name_count+=1	
	print "%2d. %s " %(file_name_count, file_name)
print "------------------------------------------------------------------------------"


for file_name in file_name_list:

	columname_line=1
	if(figure_ID==11 or figure_ID==12 or RAW_SCORE_FILE): columname_line=0 #special case since silent_file is created with GREP SCORE, no sequence line. #Mod out on March 28, 2011

	if(is_non_empty_silent_file( file_name )==False): 
		print "file_name (%s) is empty! " %(file_name )
		continue

	min_score=-100
	max_score=100
	#(min_score,max_score)=get_min_max(file_name) #TIME CONSUMING...
	if(min_score<global_min): global_min=min_score
	if(max_score>global_max): global_max=max_score

	firstlines = popen_and_readlines('head -n 3 '+ file_name, Is_master=True )
	col_name_list=firstlines[columname_line].split() 

	try:
		rmsd_column=col_name_list.index(rmsd_colname)+1
	except:
		error_exit_with_message("Cannot find score_col_index with corresponding rmsd_colname= %s in col_name_list: %s" %(rmsd_colname, list_to_string(col_name_list)) )


	if(one_point_type==False): 
		file_ID+=1
		point_type_ID=file_ID


	grep_filename_command='\'%s\'' %(file_name)
	#if(quick and (recal_rmsd==False) ): grep_filename_command= ' \'<grep "SCORE:" %s \' ' %(file_name) #Mod out on March 28, 2011
	grep_filename_command= ' \'<grep "SCORE:" %s \' ' %(file_name)

	plot_command+=" %s using %d:%d ls %d," % (grep_filename_command ,rmsd_column,score_column, point_type_ID) 

plot_command=plot_command[:-1] #Remove the last comma

#AUTO:
outfile.write('set yrange[%d:%d]\n' %(global_min-1, global_min+14))
outfile.write('set xrange[0:10]\n')

'''
# 1S72 GAAA tetraloop FARFAR weight
outfile.write('set yrange[-64:-50]\n')
outfile.write('set xrange[0:10]\n')

# 1S72 GAAA tetraloop SWA weight
outfile.write('set yrange[*:-30]\n')
outfile.write('set xrange[0:10]\n')
'''

if(figure_ID==2):
	outfile.write('set yrange[-66:-51]\n')
	outfile.write('set xrange[0:10]\n')
	FARFAR_command_before=" '/Users/sripakpa/minirosetta/10_2010_SWA_PAPER_BENCHMARK/Oct_25_SWA_1S72_RANDOMNESS_CHECK_RERUN/REGION_1_0/SCORE_start_from_region_2_0_sample.out' using 27:2 ls 7, '/Users/sripakpa/minirosetta/10_2010_SWA_PAPER_BENCHMARK/Oct_25_SWA_1S72_RANDOMNESS_CHECK_RERUN/REGION_2_1/SCORE_start_from_region_3_1_sample.out' using 27:2 ls 7, '/Users/sripakpa/minirosetta/10_2010_SWA_PAPER_BENCHMARK/Oct_25_SWA_1S72_RANDOMNESS_CHECK_RERUN/REGION_3_2/SCORE_start_from_region_4_2_sample.out' using 27:2 ls 7, '/Users/sripakpa/minirosetta/10_2010_SWA_PAPER_BENCHMARK/Oct_25_SWA_1S72_RANDOMNESS_CHECK_RERUN/REGION_4_3/SCORE_start_from_region_0_3_sample.out' using 27:2 ls 7"

################UU_UC_mismatch###############
if(figure_ID==3):
#	FARFAR_command_after=" '~/minirosetta/Nov_4_FARFAR_simple_motif_allow_bulge_mode/main_folder/2PN3_UU_UC_mismatch/recal_rmsd_cat_file.out' using 31:2 ls 6"
	FARFAR_command_after="'region_FINAL.out' using 27:2 ls 6"
##############################################

################2R8S TLR###############
if(figure_ID==4):
	outfile.write('set yrange[%d:%d]\n' %(global_min-5, global_min+84))
	outfile.write('set xrange[0:5]\n')
	FARFAR_command_after=" '~/minirosetta/10_2010/Oct_23_FARFAR_2EEW_and_MEDIUM_2R8S_DOWNLOAD/main_folder/2R8S_MEDIUM_TLR/recal_rmsd_cat_file.out' using 31:2 ls 6"
#	FARFAR_command_after="'region_FINAL.out' using 27:2 ls 6"
##############################################

################2EEW###############
if(figure_ID==5):
	outfile.write('set yrange[%d:%d]\n' %(global_min-5, global_min+54))
	outfile.write('set xrange[0:5]\n')
	FARFAR_command_after=" '~/minirosetta/10_2010/Oct_23_FARFAR_2EEW_and_MEDIUM_2R8S_DOWNLOAD/main_folder/2EEW/recal_rmsd_cat_file.out' using 31:2 ls 6"
#	FARFAR_command_after="'region_FINAL.out' using 27:2 ls 6"

if(figure_ID==6):
	outfile.write('set yrange[%d:%d]\n' %(global_min-1, global_min+5))
	outfile.write('set xrange[1:4]\n')
	FARFAR_command_after=" '~/minirosetta/10_2010/Oct_23_FARFAR_2EEW_and_MEDIUM_2R8S_DOWNLOAD/main_folder/2EEW/recal_rmsd_cat_file.out' using 31:2 ls 6"
#	FARFAR_command_after="'region_FINAL.out' using 27:2 ls 6"
##############################################

################1CSL####################
if(figure_ID==7):
	outfile.write('set yrange[%d:%d]\n' %(global_min-2, global_min+15))
	outfile.write('set xrange[0:5]\n')
	FARFAR_command_after=" '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Jan_12_500K_FARFAR/main_folder/1CSL_BULGE_SYN_CHI/FINAL_SCORE.out' using 25:2 ls 6"


################J5/Ja OLD####################
	
if(figure_ID==8):
	outfile.write('set yrange[%d:-37]\n' %(global_min-1))
	outfile.write('set xrange[0:6]\n')
	FARFAR_command_after=" '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Jan_31_FARFAR_S1_J5_Ja_2r8s_500K/main_folder/S1_2R8S_J5_5a/FINAL_SCORE.out' using 25:2 ls 6"

################2PN4####################
if(figure_ID==9):
	optimize_native_RMSD=1.492
	optimize_native_Energy=28.021

	optimize_native_filename="optimize_native_data.txt"
	if(exists(optimize_native_filename)): 
		print "optimize_native_filename (%s) already exist! ...removing" %(optimize_native_filename)
		submit_subprocess("rm %s " %(optimize_native_filename) )

	OPTIMIZE_NATIVE_FILE = open( optimize_native_filename, 'w' ) 
	OPTIMIZE_NATIVE_FILE.write('%s		%s' %(optimize_native_RMSD, optimize_native_Energy) )
	OPTIMIZE_NATIVE_FILE.close()

	FARFAR_command_after=" '%s' using 1:2 ls 6" %(optimize_native_filename) 



	outfile.write('set yrange[21:34]\n')
	outfile.write('set xrange[0:7]\n')

################J5/Ja NEW####################
if(figure_ID==10):
	outfile.write('set yrange[%d:-37]\n' %(global_min-1))
	outfile.write('set xrange[0:6]\n')
	FARFAR_command_after=" '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Mar_10_250K_FARFAR_LOOP/PART_2/main_folder/S1_2R8S_J5_5a/FINAL_SCORE.out' using 25:2 ls 6"

#################3D2V###############
if(figure_ID==11):
	outfile.write('set yrange[%d:%d]\n' %(-95.107-1,-70))
	#outfile.write('set yrange[%d:-70]\n' %(-95.107-1))
	#outfile.write('set yrange[%d:-85]\n' %(global_min-1))
	outfile.write('set xrange[0:10]\n')
	FARFAR_command_after=" '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/Mar_10_250K_FARFAR_LOOP/PART_2/main_folder/LOOP_3D2V_117_121/FINAL_SCORE.out' using 25:2 ls 31"
	
if(figure_ID==12):
	outfile.write('set yrange[%d:%d]\n' %(-95.107-1,-70))
	outfile.write('set xrange[0:10]\n')
	FARFAR_command_after=" '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/03_2011/Mar_10_250K_FARFAR_LOOP/PART_2/main_folder/LOOP_3D2V_117_121/FINAL_SCORE.out' using 25:2 ls 31"

#############A-A platform####################
		

if(figure_ID==14):
	outfile.write('set yrange[%d:%d]\n' %(-50,-34))
	outfile.write('set xrange[0:8]\n')

if(FARFAR_command_before!=""):
	plot_command=plot_command.replace("plot ", "plot %s ," %(FARFAR_command_before))

if(FARFAR_command_after!=""):
	plot_command+= ", " + FARFAR_command_after + " "

if(figure_ID==13):
	outfile.write('set yrange[-140:-100]\n') 
	outfile.write('set xrange[0:14]\n')
	plot_command="plot 'FINAL_SCORE.out' using 25:2 ls 7"


###Hack###
#if(figure_ID==12):
#	outfile.write('set yrange[%d:%d]\n' %(-95.107-1,-70))
#	outfile.write('set xrange[0:10]\n')
#	FARFAR_command_before=" '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/03_2011/Mar_10_250K_FARFAR_LOOP/PART_2/main_folder/LOOP_3D2V_117_121/FINAL_SCORE.out' using 25:2 ls 31"
#	plot_command="plot %s " %(FARFAR_command_before)
##########

outfile.write(plot_command + "\n")

outfile.write('set output\n')
outfile.write('reset\n\n')

outfile.close()

print "PLOTTING.....%s" %(figure_name) 
sys.stdout.flush()

submit_subprocess('gnuplot %s' % script_name)
submit_subprocess('gv %s &' % figure_name)

#system('open -a preview %s' %figure_name)


