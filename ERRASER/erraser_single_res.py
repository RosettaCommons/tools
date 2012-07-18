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
print 'Starting erraser_rebuild.py...'	
###########Input options####################
start_time=time.time()
input_pdb = parse_options(sys.argv, 'pdb', '')
out_prefix = parse_options(sys.argv, 'out_prefix', basename(input_pdb).replace('.pdb', '_rebuild') )
map_file = parse_options(sys.argv, 'map', '')
map_reso = parse_options(sys.argv, 'map_reso', 2.0)
verbose = parse_options(sys.argv, 'verbose', "False")
debug = parse_options(sys.argv, 'debug', "False")
new_torsional_potential= parse_options( sys.argv, "new_torsional_potential", "True" )
native_screen_RMSD= parse_options(sys.argv, "native_screen_RMSD", 2.0)
native_edensity_cutoff= parse_options(sys.argv, "native_edensity_cutoff", 0.9) 
ideal_geometry =  parse_options( sys.argv, "ideal_geometry", "True" )
include_native =  parse_options( sys.argv, "include_native", "False" )
finer_sampling = parse_options( sys.argv, "finer_sampling", "False" )
kept_temp_folder = parse_options ( sys.argv, 'kept_temp_folder', 'False' )
rebuild_res = parse_options ( sys.argv, 'rebuild_res', "" )

if input_pdb == "" : 
    error_exit("USER need to specify -pdb option")
check_path_exist( input_pdb )

if rebuild_res == "" : 
    error_exit("USER need to specify -rebuild_res option")

if map_file == "" : 
    error_exit("USER need to specify -map option")
check_path_exist( map_file )

#######File Paths#######################################################
python_file_path = os.path.split( os.path.abspath(__file__) ) [0]
SWA_rebuild_python = "%s/SWA_rebuild_erraser.py" % python_file_path
check_path_exist( SWA_rebuild_python )
minimize_python = "%s/erraser_minimize.py" % python_file_path
check_path_exist(minimize_python)
input_pdb = abspath(input_pdb)
map_file = abspath(map_file)
out_prefix = abspath(out_prefix)

if debug :
    verbose = True
    kept_temp_folder = True
#####Set temp folder#######################
base_dir = os.getcwd()
temp_dir = '%s/%s/' % (base_dir, basename(input_pdb).replace('.pdb', '_single_res_rebuild_temp') )
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

regularized_input_pdb = basename(input_pdb).replace('.pdb', '_regularized.pdb')
regularize_pdb(input_pdb, regularized_input_pdb)
[res_conversion_list, fixed_res_final, cutpoint_final, CRYST1_line] = pdb2rosetta(regularized_input_pdb, 'start.pdb')
rna_rosetta_ready_set('start.pdb', 'temp.pdb')
rebuild_res_rosetta = res_conversion_list.index(rebuild_res) + 1
cutpoint_final.sort()

print 'Starting to rebuild residue %s' % rebuild_res

command = '%s -pdb temp.pdb' % SWA_rebuild_python
if map_file != '' :
    command += ' -map '+ map_file 
    command += ' -map_reso %f' % map_reso
command += ' -ideal_geometry %s' % ideal_geometry
command += ' -include_native %s' % include_native
command += ' -slice_nearby True'
command += ' -finer_sampling %s' % finer_sampling
command += ' -rebuild_res %d' % rebuild_res_rosetta    
command += ' -native_screen_RMSD %s' % native_screen_RMSD
command += ' -native_edensity_cutoff %s' % native_edensity_cutoff
command += ' -verbose %s' % verbose
command += ' -new_torsional_potential %s ' % new_torsional_potential
command += ' -cutpoint_open '
for cutpoint in cutpoint_final :
    command += '%d ' % cutpoint
print command
subprocess_call(command, 'rebuild_%s.out' % rebuild_res, 'rebuild_%s.err' % rebuild_res)

rebuilt_pdb_merge0 = './temp_pdb_res_%d/output_pdb/S_000000_merge.pdb' % rebuild_res_rosetta
if not exists(rebuilt_pdb_merge0) :
    print "No suitable alternative structure can be sampled."
else :
    print "Residue %s is sucessfully rebuilt!" % rebuild_res

    command = minimize_python
    if map_file != '' :
        command += ' -map '+ map_file 
        command += ' -map_reso %f' % map_reso
    command += ' -ideal_geometry %s' % ideal_geometry
    command += ' -rebuild_res %d' % rebuild_res_rosetta    
    command += ' -new_torsional_potential %s ' % new_torsional_potential
    command += ' -res_slice %s ' % rebuild_res_rosetta

    final_pdb_list = []
    for i in xrange(0,5) :
        rebuilt_pdb_merge = './temp_pdb_res_%d/output_pdb/S_00000%d_merge.pdb' % (rebuild_res_rosetta, i)
        if exists(rebuilt_pdb_merge) :
            min_pdb = rebuilt_pdb_merge.replace('merge.pdb', 'min.pdb')
            min_out = './temp_pdb_res_%d/output_pdb/minimize_%d.out' % (rebuild_res_rosetta, i)
            min_err = './temp_pdb_res_%d/output_pdb/minimize_%d.err' % (rebuild_res_rosetta, i)
            minimize_command = command
            minimize_command += " -pdb %s " % rebuilt_pdb_merge
            minimize_command += " -out_pdb %s " % min_pdb
            subprocess_call(minimize_command, min_out, min_err)

            score = 0.0
            min_out_lines = open(min_out).readlines()
            for j in xrange( len(min_out_lines) - 1, -1, -1) :
                if "Total weighted score:" in min_out_lines[j] :
                    score = float(min_out_lines[j].split()[-1])
                    break
            final_pdb_list.append([score, min_pdb])
        else :
            break

final_pdb_list = sorted(final_pdb_list)
out_score = open("../scores.out" ,'w')
out_score.write("rebuilt_model score\n")
for i in xrange(0, len(final_pdb_list) ) :
    out_score.write("%d %.3f\n" % (i, final_pdb_list[i][0]) )
    rosetta2phenix_merge_back(regularized_input_pdb, final_pdb_list[i][1], "%s_%d.pdb" % (out_prefix, i) )
out_score.close()

if not verbose :
    remove('temp_pdb_res_%d' % rebuild_res_rosetta)

print '###################################'	

if not kept_temp_folder :
    os.chdir(base_dir)
    remove(temp_dir) 

total_time=time.time()-start_time
print '\n', "DONE!...Total time taken= %f seconds" % total_time
print '###################################'	

    


