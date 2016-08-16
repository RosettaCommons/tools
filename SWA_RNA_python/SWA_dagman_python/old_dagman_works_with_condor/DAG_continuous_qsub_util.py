#!/usr/bin/env python

from sys import argv, exit, stdout
import string
from glob import glob
from os import system, popen, getcwd
from os.path import basename, exists, expanduser, dirname, abspath
import os
import time
from time import sleep
###############################################################################
from SWA_dagman_python.utility.SWA_util import *
###############################################################################

SLAVE_EXE = get_PYEXE('dagman/DAG_slave.py')

SEPERATE_CLUSTERER_SLAVE_NODES = False
NUM_CLUSTERER_SLAVE_NODE = -1
print "SEPERATE_CLUSTERER_SLAVE_NODES=%s | NUM_CLUSTERER_SLAVE_NODE=%d" %(SEPERATE_CLUSTERER_SLAVE_NODES, NUM_CLUSTERER_SLAVE_NODE)
if(SEPERATE_CLUSTERER_SLAVE_NODES==True and NUM_CLUSTERER_SLAVE_NODE<=0): error_exit_with_message("SEPERATE_CLUSTERER_SLAVE_NODES==True and NUM_CLUSTERER_SLAVE_NODE<=0")

###############################################################################
def wait_until_clusterer_slave_nodes_run():

	if(SEPERATE_CLUSTERER_SLAVE_NODES==False): error_exit_with_message("SEPERATE_CLUSTERER_SLAVE_NODES==False!")

	while( True ):

		qjobs_lines = popen_and_readlines('qjobs -w | grep " RUN " ', True)
		running_clusterer_slave_nodes = []

		for n in range( NUM_CLUSTERER_SLAVE_NODE ):

			slave_dir = 'SLAVE_JOBS/%d' % n 
			slave_tag = abspath( slave_dir ).replace('/','_')

			for line in qjobs_lines:
				if( line.find( slave_tag+' ' )):
					running_clusterer_slave_nodes.append(slave_tag)
					break

		if(len(running_clusterer_slave_nodes)>=(NUM_CLUSTERER_SLAVE_NODE/2)): #At least half of them are running!
			print "At least half (%d) of the clusterer_slave_nodes are running!" %(NUM_CLUSTERER_SLAVE_NODE/2)
			for n in range(len(running_clusterer_slave_nodes)):
				print "running_clusterer_slave_nodes #%d: %s " %(n+1,  running_clusterer_slave_nodes[n] )
			return
		else:
			print "Waiting for clusterer_slave_nodes to RUN, SO_FAR %d are running" %(len(running_clusterer_slave_nodes))

		sleep(2)




