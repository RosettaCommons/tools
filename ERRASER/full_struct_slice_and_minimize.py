#!/usr/bin/env python

from os.path import abspath, exists, basename
from sys import argv
import os
import shutil
import os.path
import time
import imp

try :
    import erraser_util
except :
    file_path = os.path.split( os.path.abspath(__file__) ) [0]
    imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

print '###################################'	
print 'Starting full_struct_slice_and_minimize.py...'	
start_time=time.time()

###Load in cmdline options####
input_pdb = parse_options(argv, 'pdb', '')
map_file = parse_options(argv, 'map', '')
out_pdb = parse_options(argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_full_minimize.pdb') )
map_reso = parse_options(argv, 'map_reso', 2.0)
vary_geometry= parse_options( argv, "vary_geometry", "True" )
new_torsional_potential= parse_options( argv, "new_torsional_potential", "True" )
kept_temp_folder = parse_options ( argv, 'kept_temp_folder', 'False' )
fixed_res = parse_option_int_list ( argv, 'fixed_res' )

if input_pdb == "" : 
    error_exit_with_message("USER need to specify -pdb option")
check_path_exist(input_pdb)

if map_file != "" : 
    check_path_exist( map_file )
    map_file = abspath( map_file )

if exists(out_pdb) :
    print "Output pdb file %s exists... Remove it..." % out_pdb
    os.remove(out_pdb)
#########Exe paths#########################
python_file_path = os.path.split( os.path.abspath(__file__) ) [0]
minimize_python = "%s/erraser_minimize.py" % python_file_path
check_path_exist(minimize_python)
input_pdb = abspath(input_pdb)
out_pdb = abspath(out_pdb)
#####Set temp folder#######################
base_dir = os.getcwd()
temp_dir = '%s/%s' % (base_dir, basename(input_pdb).replace('.pdb', '_full_minimize_temp') )
if exists(temp_dir) :
    print 'Temporay directory %s exists... Remove it and create a new folder.' % temp_dir
    shutil.rmtree(temp_dir)
    os.mkdir(temp_dir)
else :
    print 'Create temporary directory %s...' % temp_dir
    os.mkdir(temp_dir) 
###########################################
os.chdir(temp_dir)
shutil.copy( input_pdb, 'before_min.pdb' ) 
total_res = get_total_res(input_pdb)
n_chunk = int(total_res / 100.0 + 0.5)
cmdline_common = minimize_python
cmdline_common += ' -pdb before_min.pdb'
cmdline_common += ' -out_pdb after_min.pdb'
if map_file != '' :
    cmdline_common += ' -map %s' % map_file
    cmdline_common += ' -map_reso %s' % map_reso
cmdline_common += ' -vary_geometry %s' % vary_geometry
cmdline_common += ' -new_torsional_potential %s' % new_torsional_potential
cmdline_common += ' -fixed_res'
for res in fixed_res :
    cmdline_common += ' %d' % res
############################################
if n_chunk <= 1 :
    print "Input pdb < 150 residus, no slicing is required."
    print "Start minimizing the full pdb..."
    cmdline = cmdline_common
    cmdline += " > full_minimize_temp.out 2> full_minimize_temp.err"
    subprocess_call(cmdline)
    subprocess_call('mv after_min.pdb %s' % out_pdb)
else :
    print "Input pdb >= 150 residus, slice into %s chunks and minimize each one sequentially." % n_chunk
    res_slice_list = pdb_slice_into_chunks(input_pdb, n_chunk)
    current_chunk = 0
    for res_slice in res_slice_list :
        current_chunk += 1
        cmdline = cmdline_common
        cmdline += ' -res_slice'
        for res in res_slice :
            cmdline += ' %d' % res
        cmdline += " > full_minimize_temp_%s.out 2> full_minimize_temp_%s.err" % (current_chunk, current_chunk)
        print "Start minimizing chunk %s..." % current_chunk
        print "Chunk %s residues: %s" % (current_chunk, res_slice)
        subprocess_call(cmdline)
        shutil.copy('after_min.pdb', 'before_min.pdb')
        print "Minimization for chunk %s ends sucessfully." % current_chunk
    subprocess_call('mv after_min.pdb %s' % out_pdb)

if not kept_temp_folder :
    os.chdir(base_dir)
    shutil.rmtree(temp_dir) 

total_time=time.time()-start_time

print '\n', "DONE!...Total time taken= %f seconds" % total_time
print '###################################'	
