#!/usr/bin/env python
import os.path
import imp

try :
    import erraser_util
except :
    file_path = os.path.split( os.path.abspath(__file__) ) [0]
    imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

print '###################################'	
print 'Starting seq_rebuild.py...'	
###########Input options####################
start_time=time.time()
input_pdb = parse_options(sys.argv, 'pdb', '')
out_pdb = parse_options(sys.argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_rebuild.pdb') )
map_file = parse_options(sys.argv, 'map', '')
map_reso = parse_options(sys.argv, 'map_reso', 2.0)
verbose = parse_options(sys.argv, 'verbose', "False")
debug = parse_options(sys.argv, 'debug', "False")
new_torsional_potential= parse_options( sys.argv, "new_torsional_potential", "True" )
native_screen_RMSD= parse_options(sys.argv, "native_screen_RMSD", 2.0)
native_edensity_cutoff= parse_options(sys.argv, "native_edensity_cutoff", 0.9) 
ideal_geometry =  parse_options( sys.argv, "ideal_geometry", "True" )
include_native =  parse_options( sys.argv, "include_native", "False" )
slice_nearby =  parse_options( sys.argv, "slice_nearby", "True" )
finer_sampling = parse_options( sys.argv, "finer_sampling", "False" )
kept_temp_folder = parse_options ( sys.argv, 'kept_temp_folder', 'False' )
rebuild_res = parse_option_int_list ( sys.argv, 'rebuild_res' )
cutpoint_open = parse_option_int_list ( sys.argv, 'cutpoint_open' )

if input_pdb == "" : 
    error_exit("USER need to specify -pdb option")
check_path_exist( input_pdb )

if len(rebuild_res) == 0 : 
     error_exit("USER need to specify -rebuild_res option")

if map_file != "" : 
    check_path_exist( map_file )
    map_file = abspath( map_file )

check_path_exist( input_pdb )
if exists(out_pdb) :
    print "Output pdb file %s exists... Remove it..." % out_pdb
    remove(out_pdb)

if debug :
    verbose = True
    kept_temp_folder = True
###########Exe file paths###########################
python_file_path = os.path.split( os.path.abspath(__file__) ) [0]
SWA_rebuild_python = "%s/SWA_rebuild_erraser.py" % python_file_path
check_path_exist( SWA_rebuild_python )
input_pdb = abspath(input_pdb)
out_pdb = abspath(out_pdb)
#####Set temp folder#######################
base_dir = os.getcwd()
temp_dir = '%s/%s/' % (base_dir, basename(input_pdb).replace('.pdb', '_seq_rebuild_temp') )
if exists(temp_dir) :
    print 'Temporay directory %s exists... Remove it and create a new folder.' % temp_dir
    remove(temp_dir)
    os.mkdir(temp_dir)
else :
    print 'Create temporary directory %s...' % temp_dir
    os.mkdir(temp_dir) 

print '###################################'	
#####################################################
os.chdir(temp_dir)

copy(input_pdb, 'temp.pdb')

total_res = get_total_res(input_pdb)
sucessful_res = []
failed_res = []
for res in rebuild_res:
    print 'Starting to rebuild residue %s' % res

    command = '%s -pdb temp.pdb' % SWA_rebuild_python
    if map_file != '' :
        command += ' -map '+ map_file 
        command += ' -map_reso %f' % map_reso
    command += ' -ideal_geometry %s' % ideal_geometry
    command += ' -include_native %s' % include_native
    command += ' -slice_nearby %s' % slice_nearby
    command += ' -finer_sampling %s' % finer_sampling
    command += ' -rebuild_res %d' % res    
    command += ' -native_screen_RMSD %s' % native_screen_RMSD
    command += ' -native_edensity_cutoff %s' % native_edensity_cutoff
    command += ' -verbose %s' % verbose
    command += ' -new_torsional_potential %s ' % new_torsional_potential
    command += ' -cutpoint_open '
    for cutpoint in cutpoint_open :
        command += '%d ' % cutpoint
    print command
    subprocess_call(command, 'seq_rebuild_temp_%s.out' % res, 'seq_rebuild_temp_%s.err' % res)

    if exists('./temp_pdb_res_%d/output_pdb.txt' % res) :
        print 'Job completed for residue %s' % res

        rebuilt_pdb_merge = './temp_pdb_res_%d/output_pdb/S_000000_merge.pdb' % res
        rebuilt_pdb_orig = './temp_pdb_res_%d/output_pdb/S_000000.pdb' % res

        if exists(rebuilt_pdb_merge) :
            copy(rebuilt_pdb_merge, 'temp.pdb')
            sucessful_res.append(res)
            print "Residue %d is sucessfully rebuilt!" % res
        elif exists(rebuilt_pdb_orig) :
            copy(rebuilt_pdb_orig, 'temp.pdb')
            sucessful_res.append(res)
            print "Residue %d is sucessfully rebuilt!" % res
        else :
            failed_res.append(res)
            print "No suitable alternative structure can be sampled."
            print "Residue %d is not rebuilt!" % res
    else :
        error_exit("ERROR in rebuilding!!!")
    if not verbose :
        remove('temp_pdb_res_%d' % res)

    print '###################################'	
    
copy('temp.pdb', out_pdb)

if not kept_temp_folder :
    os.chdir(base_dir)
    remove(temp_dir) 

print "All rebuilding moves completed sucessfully!!!!"
print 'sucessful_res: %s' % sucessful_res
print 'failed_res: %s' % failed_res

total_time=time.time()-start_time
print '\n', "DONE!...Total time taken= %f seconds" % total_time
print '###################################'	

    


