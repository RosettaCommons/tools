#!/usr/bin/env python
"""
This script removes neighboring points from free or working sets

Usage: layerline_pruner.py -r layer_lines1 -t layer_lines2 -d delta

layer_lines1 - reference point,
layers_lines2 - o be prunned set
delta -running distance in reciprocal space 
(for points seperated at 0.0025, delta=0.0035 should work)  
"""
__author__ = "Wojtek Potrzebowski"
__credits__ = ["Wojtek Potrzebowski and Ingemar Andre"]
__maintainer__ = "Wojtek Potrzebowski"
__email__ = "Wojciech.Potrzebowski@biochemistry.lu.se"

import optparse
import os
import sys 
import numpy as np

def read_diff_file(filename):
	f = open(filename)
	layer_lines_start = {}
	layer_lines_end = {}
	R_l = {}
	I_l = {}

	lines = f.readlines()
	last_index=0
	lindex= 0
	layer_index = 0
	layer_lines_start[0] = float(lines[0].split(" ")[1])
	R_l[0] = [layer_lines_start[0]]
	I_l[0] = [float(lines[0].split(" ")[0])]

	for line in lines:
		line_arr = line.split(" ")
		layer_index = int(line_arr[2].strip())
		if layer_index!=last_index:
			layer_lines_end[layer_index-1] = float(lines[lindex-1].split(" ")[1])
			layer_lines_start[layer_index] = float(line_arr[1])
			R_l[layer_index] = [layer_lines_start[layer_index]]
			I_l[layer_index] = [float(line_arr[0])]
			last_index = layer_index
		else:
			R_l[layer_index].append(float(line_arr[1]))
			I_l[layer_index].append(float(line_arr[0]))
		lindex+=1
	return R_l, I_l

def remove_close_points(WR,FR,FI, delta):
	"""
	Removve points in R space which are closer than delta in work set
	"""
	trimed_FR = {}
	trimed_FI = {}
	removed_FR = {}
        removed_FI = {}
	for lline in FR.keys():
		lindex = 0
		trimed_FR[lline]=[]
		trimed_FI[lline]=[]
		removed_FR[lline]=[]
                removed_FI[lline]=[]
		for fR in FR[lline]:
			remove = 0
			for wR in WR[lline]:
				if abs(fR-wR)<=delta:
					#if lindex+1 %2 == 0:
					#	print "Removed", fR, wR
					remove=1
			if remove==1:
				fI = FI[lline][lindex]
                                removed_FR[lline].append(fR)
                                removed_FI[lline].append(fI)
			if remove==0:
				fI = FI[lline][lindex]
				trimed_FR[lline].append(fR)
				trimed_FI[lline].append(fI)
			lindex+=1
	return trimed_FR, trimed_FI, removed_FR, removed_FI

def write_to_file_removed(free_set, trimed_FR, trimed_FI):
        out_file = open("rem_"+free_set,"w")
        for lline in trimed_FR.keys():
                lindex = 0
                for fR in trimed_FR[lline]:
                        out_file.write(str(trimed_FI[lline][lindex])+" "+str(fR)+" "+str(lline)+"\n")
                        lindex+=1
        out_file.close()                
                        

def write_to_file(free_set, trimed_FR, trimed_FI):
	out_file = open("trm_"+free_set,"w")
	for lline in trimed_FR.keys():
                lindex = 0
                for fR in trimed_FR[lline]:
			out_file.write(str(trimed_FI[lline][lindex])+" "+str(fR)+" "+str(lline)+"\n")
			lindex+=1
	out_file.close()		
			

if __name__=="__main__":
	doc = """
	    Prunes points from target layer lines with the respect to reference layer lines.
	    Program outputs prunned layer lines (trm prefix) and removed points (rem). 
	    Usage: python layer_line_pruner.py --help
	"""
	print doc
	usage = "usage: %prog [options] args"
	option_parser_class = optparse.OptionParser
	parser = option_parser_class( usage = usage, version='0.1' )

    	parser.add_option("-r", "--reference", dest="reference", 
                      help="Reference layer lines [OBLIGATORY]")
   	parser.add_option("-t", "--target", dest="target",
                      help="To be prunned layer lines [OBLIGATORY]")
    	parser.add_option("-d", "--delta", dest="delta", default =0.0035,
                      type = 'float',
                      help="Prunning distance, e.g. 0.0035 for points sepoarated at 0.0025  [OBLIGATORY]")
    	options, args = parser.parse_args()
	
	work_set = options.reference
	free_set = options.target
	delta = options.delta
	WR, WI = read_diff_file(work_set)
	FR, FI = read_diff_file(free_set)
	trimed_FR, trimed_FI, removed_FR, removed_FI = remove_close_points(WR,FR,FI, delta)
	write_to_file(free_set, trimed_FR, trimed_FI)
	write_to_file_removed(free_set, removed_FR, removed_FI)	
#Last line
