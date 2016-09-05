#!/usr/bin/env python

from os import system,popen
from os.path import exists,dirname,basename,expanduser
from sys import exit, argv
import string
import os
from time import sleep

filename=argv[1]
type_col=1 #0
num_col=0 #26 


infile=open(filename,"r")


desired_seq_string=argv[2]
print "desired_seq_string: %s" %(desired_seq_string)


desired_seq_list=desired_seq_string.split('-')
print "desired_seq_list: ", desired_seq_list


pyrimidine=False

seq_list=[]

for line in infile:
#	print line

	line_list=line.split()
	if(len(line_list)<=type_col): continue
	if(len(line_list)<=num_col): continue


#	print line_list
	res_type=line_list[type_col]
	res_num=int(line_list[num_col])

	seq=(res_type, res_num)
#	print seq
	seq_list.append(seq)

#for i in range(len(a))

for n in range(len(seq_list)):
	if(n>len(seq_list)-len(desired_seq_list)): continue

	found_match=True

	for i in range(len(desired_seq_list)):
#		print seq_list[n+i][0], desired_seq_list[i], (seq_list[n+i][0]!=desired_seq_list[i])
		if(desired_seq_list[i]=='x'): continue

		if(seq_list[n+i][0]!=desired_seq_list[i]): 
			found_match=False
			break

	if(found_match==True):
		check_string=''
		for i in range(len(desired_seq_list)):			
			check_string+="(%s,%d)" % (seq_list[n+i][0],seq_list[n+i][1]) 
			
		print check_string


