#!/usr/bin/env python
import os.path
import imp

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

#####################################################
print '###################################'
print 'Starting RNA_rosetta_ready_set.py...'
start_time=time.time()
#####Input options#####################################
input_pdb = parse_options(sys.argv, 'pdb', '')
out_pdb = parse_options(sys.argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_ready_set.pdb') )

if input_pdb =="" :
    error_exit("Error: USER need to specify -pdb option")
check_path_exist(input_pdb)

if exists( out_pdb ) :
    print 'Output pdb %s already exists... Remove it.' % out_pdb
    remove(out_pdb)
#######File Paths#######################################################
base_dir = os.getcwd()
temp_file = '%s/%s_temp_rs.pdb' % (base_dir, basename(input_pdb).replace('.pdb', '') )
input_pdb = abspath(input_pdb)
out_pdb = abspath(out_pdb)

pdb2rosetta(input_pdb, temp_file)
rna_rosetta_ready_set(temp_file, out_pdb)
remove(temp_file)

total_time=time.time()-start_time
print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
print '###################################'

