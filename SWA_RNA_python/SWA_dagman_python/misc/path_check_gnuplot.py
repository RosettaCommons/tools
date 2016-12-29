#!/usr/bin/env python

from os import system,popen,path
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
from time import sleep
import os

assert( len(argv)>1)

num_elements=int(argv[1])

figure_name=os.path.abspath(os.curdir)
figure_name=figure_name.replace('/','_')
figure_name=figure_name + '_check_path.eps'

print "figure_name= %s" %figure_name


script_name='check_path_script.txt'

system( 'rm %s' % script_name)  #create parent directory as needed
system( 'rm %s' % figure_name)  #create parent directory as needed

outfile = open( script_name, 'w' ) 

outfile.write('\nset terminal postscript enhanced color\n')
outfile.write('set output "%s"\n' % figure_name)
outfile.write('set style line 7 lt 3 lw 4\n')
outfile.write('set style line 5 lt 1 lw 4\n')
outfile.write('set style line 8 lt 3 lw 4\n')
outfile.write('set style line 9 lt 1 lw 4\n')

outfile.write('set style line 1 linetype "solid" lw  5\n')
outfile.write('set border lt rgb "yellow" lw 1\n')

outfile.write('set nokey\n')
outfile.write('set size 1.0, 1.0\n')
outfile.write('set pointsize 0.4\n')

outfile.write('unset xtics\n')
outfile.write('unset ytics\n')
#outfile.write('set xtics (1.0)\n')
outfile.write('set format x ""\n')
outfile.write('set format y ""\n')
outfile.write('set grid  lw 5\n') 

outfile.write('set yrange[-3:3]\n')
outfile.write('set xrange[-3:3]\n')
outfile.write('set multiplot\n')


screen_width=1.0
screen_height=1.0;
rmsd_column=22
score_column=2

i_missing=0
j_missing=0

if(len(argv)>3):
	i_missing=int(argv[2])
	j_missing=int(argv[3])

i_start=0
j_start=0	

if(len(argv)>5):
	i_start=int(argv[4])
	j_start=int(argv[5])


plot_width=(screen_width)/(num_elements-j_missing)
plot_height=(screen_height*0.95)*(num_elements)/((num_elements-i_missing)*(num_elements))

main_job_folder=""
extension=""
if(exists('./DAG_JOB/')):
	main_job_folder='./DAG_JOB/'
	extension='dag_job'
elif(exists('./CONDOR/')): 
	main_job_folder='./CONDOR/'
	extension='condor'
else:
	print "cannot find main_job_folder!"
	assert(False)

if( exists("%s/SAMPLER/" %(main_job_folder) ) ):
	main_job_folder="%s/SAMPLER/" %(main_job_folder)

for i in range( num_elements ) : #column
	for j in range( num_elements ) : #row

		x_coord=((plot_width)*(j-j_start))

		col= ((i)) % num_elements

		if(col==0): col=num_elements

		y_coord=screen_height*0.975-((plot_height)*(col-i_start))

		folder_name="%s/REGION_%d_%d" % (main_job_folder, i,j)
		file_name1=folder_name + "_START_FROM_REGION_%d_%d.%s" % ((i+1) % num_elements , (j) % num_elements,extension)
		file_name2=folder_name + "_START_FROM_REGION_%d_%d.%s" % ((i) % num_elements , (j-1) % num_elements,extension)
		file_name3=folder_name + "_START_FROM_REGION_%d_%d.%s" % ((i+2) % num_elements , (j) % num_elements,extension)
		file_name4=folder_name + "_START_FROM_REGION_%d_%d.%s" % ((i) % num_elements , (j-2) % num_elements,extension)
		file_name5=folder_name + "_START_FROM_REGION_0_%d_AND_%d_0.%s" % ((j) % num_elements , (i+1) % num_elements,extension)
		file_name6=folder_name + "_START_FROM_REGION_0_%d_AND_%d_0.%s" % ((j-1) % num_elements , (i) % num_elements,extension)

		outfile.write('\nunset label\n')

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

		"""
		#HACKY
		file_name1=""
		file_name2=""
		file_name3=""
		file_name4=""
		file_name5=folder_name + "_START_FROM_REGION_0_%d_AND_%d_0.%s" % ((j) % num_elements , (i+2) % num_elements,extension)
		file_name6=folder_name + "_START_FROM_REGION_0_%d_AND_%d_0.%s" % ((j-2) % num_elements , (i) % num_elements,extension)
		"""

#		file_name_cluster="region_%d_%d_sample.cluster.out" % (i,j)

		if ( ( exists( file_name1 ) or exists( file_name2 ) or exists( file_name3 ) or exists( file_name4 ) or exists( file_name5 ) or exists( file_name6 ) )==False ): continue

		print "file_name1=%s " %(file_name1) 
		print "file_name2=%s " %(file_name2) 
		print "file_name3=%s " %(file_name3) 
		print "file_name4=%s " %(file_name4) 
		print "file_name5=%s " %(file_name5) 
		print "file_name6=%s " %(file_name6) 



		print "i= %d j= %d "  % (i, j)
		outfile.write('\nunset label\n')
		outfile.write('set label "%d %d" at graph 0.4, graph -0.08 \n' % (i, j) )
		outfile.write('set size %8.4f , %8.4f\n' % (plot_width, plot_height) ) 
		outfile.write('set origin %8.4f , %8.4f\n' % (x_coord, y_coord) )


		plot_command='plot '
		if(exists( file_name1 )):
			plot_command+=' 999999*x+999999 ls 5,'
		if(exists( file_name2 )):
			plot_command+=' -1 ls 9,'  
		if(exists( file_name3 )):
			plot_command+=' 999999*x-999999 ls 7,' 
		if(exists( file_name4 )):
			plot_command+=' 1 ls 8,'  
		if(exists( file_name5 )):
			plot_command+=' 10*x ls 2,'
		if(exists( file_name6 )):
			plot_command+=' 0 ls 1,'  


#		if(exists( file_name_cluster )):
#			plot_command+=' "1" using %d:%d ls 9,' % (file_name_cluster ,rmsd_column,score_column) 
		
		plot_command=plot_command[:-1] #Remove the last comma
		print "plot_command= %s" %(plot_command)

		outfile.write(plot_command + "\n")

		
outfile.write('\nunset multiplot\n')
outfile.write('set output\n')
outfile.write('reset\n\n')

outfile.close()

system('gnuplot %s' % script_name)
system('gv %s' % figure_name)


