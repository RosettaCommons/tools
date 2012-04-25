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
print 'Starting erraser_minimize.py...'	
start_time=time.time()

#######Load in cmdline options#####################
input_pdb = parse_options(sys.argv, 'pdb', '')
map_file = parse_options(sys.argv, 'map', '')
out_pdb = parse_options(sys.argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_minimize.pdb') )
map_reso = parse_options(sys.argv, 'map_reso', 2.0)
vary_geometry= parse_options( sys.argv, "vary_geometry", "True" )
constrain_phosphate = parse_options( sys.argv, "constrain_phosphate", "False" )
new_torsional_potential= parse_options( sys.argv, "new_torsional_potential", "True" )
fixed_res = parse_option_int_list ( sys.argv, 'fixed_res' )
res_slice = parse_option_int_list ( sys.argv, 'res_slice' )

if input_pdb == "" : 
    error_exit("USER need to specify -pdb option")
check_path_exist(input_pdb)

if map_file != "" : 
    check_path_exist( map_file )
    map_file = abspath( map_file )
#######Folders and files paths###########################
database_folder = rosetta_database() 
rna_minimize_exe = rosetta_bin("erraser_minimizer.linuxgccrelease")
temp_rs = input_pdb.replace('.pdb', '_temp_rs.pdb')
temp_rs_min = input_pdb.replace('.pdb', '_temp_rs_min.pdb')

if exists(temp_rs) :
    print "Temporary file %s exists... Remove it..." % temp_rs
    remove(temp_rs) 
if exists(temp_rs_min) :
    print "Temporary file %s exists... Remove it..." % temp_rs_min
    remove(temp_rs_min)
if exists(out_pdb) :
    print "Output pdb file %s exists... Remove it..." % out_pdb
    remove(out_pdb)

####slicing into smaller pdbs if given in the option####
fixed_res_final = []
res_sliced_all = []
if len(res_slice) != 0 :
    [res_sliced_all, patched_res] = pdb_slice_with_patching( input_pdb, temp_rs, res_slice )
    fixed_res_final += patched_res
    for res in fixed_res :
        if res in res_sliced_all :
            fixed_res_final.append( res_sliced_all.index(res) + 1 )
else :
    copy(input_pdb, temp_rs)
    fixed_res_final = fixed_res
    
####submit rosetta cmdline##############
command = rna_minimize_exe 
command += " -database %s " % database_folder
command += " -native %s " % temp_rs
command += " -out_pdb %s " % temp_rs_min
if map_file != '' :
    command += " -score::weights rna/rna_hires_elec_dens "
else :
    command += " -score::weights rna/rna_loop_hires_04092010"

if new_torsional_potential :
    command += " -score:rna_torsion_potential RNA09_based_2012_new "

command += " -vary_geometry %s " % str(vary_geometry).lower()
command += " -constrain_P %s " % str(constrain_phosphate).lower()

if len(fixed_res_final) != 0 :
    command += ' -fixed_res '
    for i in fixed_res_final :
        command += '%d ' % i

if map_file != '' :
    command += " -edensity:mapfile %s " % abspath(map_file)
    command += " -edensity:mapreso %s " % map_reso
    command += " -edensity:realign no "
print "cmdline: %s" % command
print "#######Submit the Rosetta Command Line###############"
subprocess_call(command)
print "#####################################################"

####Merge final result back to pdb####
if len(res_slice) != 0 :
    sliced2orig_merge_back( input_pdb, temp_rs_min, out_pdb, res_sliced_all )
    remove( temp_rs_min )
else :
    shutil.move(temp_rs_min, out_pdb)

remove( temp_rs )
#########################################

total_time=time.time()-start_time

print '\n', "DONE!...Total time taken= %f seconds" % total_time
print '###################################'	

