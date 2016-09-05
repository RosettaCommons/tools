#!/usr/bin/env python

from os import system,popen,path
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os


from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.parser.SWA_parse_options import parse_options
######################################################################



assert( len(argv)>1)




trace_pathway_folder = parse_options( argv, "trace_pathway_folder", "" )

ignore_region_list = parse_options( argv, "ignore_region_list", [""] )

rmsd_colname=parse_options(argv, "rmsd_colname", "O_loop_rmsd")

#rmsd_column=  parse_options( argv, "rmsd_column", 27 )


print "left over argv=", argv

print "rmsd_colname= %s " %(rmsd_colname)

num_elements=int(argv[1])
quick=(len(argv)==3)

i_missing=0
j_missing=0

if(len(argv)>3):
	i_missing=int(argv[2])
	j_missing=int(argv[3])
	quick=(len(argv)==5)

i_start=0
j_start=0	

if(len(argv)>5):
	i_start=int(argv[4])
	j_start=int(argv[5])
	quick=(len(argv)==7)


#screen_width=0.2*(num_elements-j_missing)
#screen_height=0.2*(num_elements-i_missing)

#screen_width=1.4
#screen_height=1.3

screen_width=1.0
screen_height=1.0



score_column=2




plot_width=(screen_width)/(num_elements-j_missing)
plot_height=(screen_height*0.95)*(num_elements)/((num_elements-i_missing)*(num_elements))

print "plot screen_width %f screen_heigth %f " %(screen_width, screen_height)
print "plot_width %f plot_height %f " %(plot_width, plot_height)

#figure_name='score_vs_rmsd.eps'

figure_name=os.path.abspath(os.curdir)
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
outfile.write('set style line 7 lt 3 lw 2\n')
outfile.write('set style line 5 lt 1 lw 2\n')

outfile.write('set style line 1 linetype "solid" lw  5\n')
outfile.write('set border lt rgb "blue" lw 1\n')

outfile.write('set nokey\n')
#outfile.write('set size 1.0, 1.0\n')
outfile.write('set size %f,%f\n' %(screen_width*1.1,screen_height*1.0) )  
outfile.write('set pointsize 0.4\n')

outfile.write('unset xtics\n')
#outfile.write('unset ytics\n')
outfile.write('set ytics font "Helvetica,10"\n')
outfile.write('set xtics (1.0)\n')
outfile.write('set format x ""\n')

#outfile.write('set format y ""\n')
outfile.write('set grid  lw 5\n') 
outfile.write('set grid front\n')


#outfile.write('set yrange[*:*]\n')
outfile.write('set xrange[0:*]\n')
#outfile.write('set xrange[0:3]\n')
#outfile.write('set multiplot\n')

#set multiplot title "Demo of placing multiple plots (2D and 3D)\nwith explicit alignment of plot borders"
outfile.write('set multiplot title \'%s\' \n' %(figure_name_inside))



    
for i in range( num_elements ) : #column
	for j in range( num_elements ) : #row

		###########################################################
		ignore_this_region=False
		if(ignore_region_list!=[""]):
			for nn in range(len(ignore_region_list)):
				ignore_region=ignore_region_list[nn].split('-')
				if(len(ignore_region)!=2): error_exit_with_message("len ignore_region_list[%d]!=2, ignore_region_list[%d]=" %(nn,nn,ignore_region_list[nn]) )
				if( i==int(ignore_region[0]) and j==int(ignore_region[1]) ): 
					ignore_this_region=True
					break

		if(ignore_this_region==True): continue
		###########################################################
			

		x_coord=((plot_width)*(j-j_start))

		col= ((i)) % num_elements

		if(col==0): col=num_elements

		y_coord=screen_height*0.975-((plot_height)*(col-i_start))


		folder_name="REGION_%d_%d" % (i,j)
		file_name1=folder_name + "/" + "start_from_region_%d_%d_sample_filtered.out" % ((i+1) % num_elements , (j) % num_elements)
		file_name2=folder_name + "/" + "start_from_region_%d_%d_sample_filtered.out" % ((i) % num_elements , (j-1) % num_elements)
		file_name3=folder_name + "/" + "start_from_region_%d_%d_sample_filtered.out" % ((i+2) % num_elements , (j) % num_elements)
		file_name4=folder_name + "/" + "start_from_region_%d_%d_sample_filtered.out" % ((i) % num_elements , (j-2) % num_elements)
		file_name5=folder_name + "/" + "start_from_region_0_%d_AND_%d_0_sample_filtered.out" % ((j) % num_elements , (i+1) % num_elements)
		file_name6=folder_name + "/" + "start_from_region_0_%d_AND_%d_0_sample_filtered.out" % ((j-1) % num_elements , (i) % num_elements)


		if(quick):
			file_name1=""
			file_name2=""
			file_name3=""
			file_name4=""
			file_name5=""
			file_name6=""

#.low40000.out
#_filtered.out"

		file_name_cluster="region_%d_%d_sample.cluster.out" % (i,j)

		trace_pathway_filename=""
		if(trace_pathway_folder!=""): trace_pathway_filename=	"%s/SCORE_%s" %(trace_pathway_folder,file_name_cluster)

#for output purposes
		is_non_empty_silent_file(file_name1, verbose=True)
		is_non_empty_silent_file(file_name2, verbose=True)
		is_non_empty_silent_file(file_name3, verbose=True)
		is_non_empty_silent_file(file_name4, verbose=True)
		is_non_empty_silent_file(file_name5, verbose=True)
		is_non_empty_silent_file(file_name6, verbose=True)
		is_non_empty_silent_file(file_name_cluster, verbose=True)
##############

		if( ( is_non_empty_silent_file(file_name1) or \
				  is_non_empty_silent_file(file_name2) or \
				  is_non_empty_silent_file(file_name3) or \
				  is_non_empty_silent_file(file_name4) or \
				  is_non_empty_silent_file(file_name5) or \
				  is_non_empty_silent_file(file_name6) or \
				  is_non_empty_silent_file(file_name_cluster) ) == False): continue

		#consistency_check
		if(exists(file_name_cluster)):
			if(is_non_empty_silent_file(file_name_cluster)==False): 
				error_exit_with_message("One of the input_silent_file of the cluster_silent_file %s is non_empty. The cluster_silent_file exists but is empty!!" %(file_name_cluster))



#		if(quick==False):
#			print "file_name1=%s " %(file_name1) 
#			print "file_name2=%s " %(file_name2) 
#			print "file_name3=%s " %(file_name3) 
#			print "file_name4=%s " %(file_name4) 
#			print "file_name5=%s " %(file_name5) 
#			print "file_name6=%s " %(file_name6) 


#		print "j= %d i= %d " % (j, i),
#		print "  x= %f y= %f "  % (x_coord, y_coord)


		outfile.write('\nunset label\n')

#		act_width=plot_width*(1+0.4*(num_elements/8.0) )
#		act_height=plot_height*(0.5+0.6*(num_elements/8.0) ) 
#		outfile.write('set size %8.4f , %8.4f\n' % (act_width , act_height ) )
#		outfile.write('set origin %8.4f , %8.4f\n' % (x_coord, y_coord) )

		width_margin=0.20
		height_margin=0.20

		lmargin=x_coord+(width_margin*plot_width)
		rmargin=x_coord+((1-width_margin)*plot_width)

		bmargin=y_coord+((height_margin)*plot_height)
		tmargin=y_coord+((1-height_margin)*plot_height)

		outfile.write('set lmargin at screen %8.4f\n' %(lmargin) )
		outfile.write('set rmargin at screen %8.4f\n' %(rmargin) )
		outfile.write('set bmargin at screen %8.4f\n' %(bmargin) )
		outfile.write('set tmargin at screen %8.4f\n' %(tmargin) )
		outfile.write('set label "%d %d" font "Helvetica,10" at graph 0.30, graph -0.15  \n' % (i, j) )



		global_min=10000000000000000
		global_max=-10000000000000000

		plot_command='plot '
		if(is_non_empty_silent_file( file_name1 )):
			(min_score,max_score)=get_min_max(file_name1)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 5,' % (file_name1 ,get_silent_file_col_index(rmsd_colname, file_name1)+1 ,score_column) 
		if(is_non_empty_silent_file( file_name2 )):
			(min_score,max_score)=get_min_max(file_name2)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 6,' % (file_name2 ,get_silent_file_col_index(rmsd_colname, file_name2)+1,score_column) 
		if(is_non_empty_silent_file( file_name3 )):
			(min_score,max_score)=get_min_max(file_name3)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 7,' % (file_name3 ,get_silent_file_col_index(rmsd_colname, file_name3)+1,score_column) 
		if(is_non_empty_silent_file( file_name4 )):
			(min_score,max_score)=get_min_max(file_name4)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 8,' % (file_name4 ,get_silent_file_col_index(rmsd_colname, file_name4)+1,score_column) 
		if(is_non_empty_silent_file( file_name5 )):
			(min_score,max_score)=get_min_max(file_name5)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 21,' % (file_name5 ,get_silent_file_col_index(rmsd_colname, file_name5)+1,score_column) 
		if(is_non_empty_silent_file( file_name6 )):
			(min_score,max_score)=get_min_max(file_name6)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 2,' % (file_name6 ,get_silent_file_col_index(rmsd_colname, file_name6)+1,score_column) 

		if(is_non_empty_silent_file( file_name_cluster )):
			(min_score,max_score)=get_min_max(file_name_cluster)
#			print "min= %f, max=%f" %(min_score,max_score)
			if(global_min>min_score): global_min=min_score
			if(global_max<max_score): global_max=max_score
			plot_command+=' "%s" using %d:%d ls 9,' % (file_name_cluster ,get_silent_file_col_index(rmsd_colname, file_name_cluster)+1,score_column) 

		if(trace_pathway_filename!=""):
			if(exists(trace_pathway_filename)):
				if(is_non_empty_silent_file(trace_pathway_filename)==False):
					error_exit_with_message("is_non_empty_silent_file(trace_pathway_filename)==False, trace_pathway_filename=%s" %(trace_pathway_filename) )
				outfile.write('set style line 1 lt 4 lw 20\n')
				plot_command+=' "%s" using %d:%d ls 1,' % (trace_pathway_filename ,get_silent_file_col_index(rmsd_colname, trace_pathway_filename)+1,score_column) 

		
#		print "global_min= %f, global_max=%f" %(global_min, global_max)
		if(global_max>global_min+30): global_max=global_min+30
			
		outfile.write('set ytics (%d,%d)\n' %(global_min, global_max-1))		
		outfile.write('set ytics offset (%8.4f,%8.4f)\n' %(0, 0.8))		
		outfile.write('set yrange[%d:%d]\n' %(global_min-1, global_max+1))

		if(is_non_empty_silent_file( file_name_cluster )):		
			outfile.write('set label "%d" font "Helvetica,10" at graph 0.30, graph +1.15 \n' %(count_struct(file_name_cluster)))

		plot_command=plot_command[:-1] #Remove the last comma
#		print "plot_command= %s" %(plot_command)

		outfile.write(plot_command + "\n")


outfile.write('\nunset multiplot\n')
outfile.write('set output\n')
outfile.write('reset\n\n')

outfile.close()

system('gnuplot %s' % script_name)
#system('gv %s' % figure_name)

system('open -a preview %s' %figure_name)


