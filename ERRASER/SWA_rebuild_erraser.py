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
print 'Starting SWA_rebuild_erraser.py...'	
start_time=time.time()
#############USER INPUT OPTIONS###########################################
input_pdb =  parse_options( sys.argv, "pdb", "" ) 
rebuild_res = parse_options( sys.argv, "rebuild_res", 0 )
map_file = parse_options( sys.argv, "map", "")
map_reso = parse_options (sys.argv, 'map_reso', 2.0)
verbose= parse_options( sys.argv, "verbose", "False" )
native_screen_RMSD= parse_options(sys.argv, "native_screen_RMSD", 2.0)
native_edensity_cutoff= parse_options(sys.argv, "native_edensity_cutoff", 0.9) 
cluster_RMSD = parse_options( sys.argv, "cluster_RMSD", 0.3 )
ideal_geometry =  parse_options( sys.argv, "ideal_geometry", "True" )
include_native =  parse_options( sys.argv, "include_native", "False" )
slice_nearby =  parse_options( sys.argv, "slice_nearby", "True" )
finer_sampling = parse_options( sys.argv, "finer_sampling", "False" )
is_append = parse_options( sys.argv, "is_append", "True" )
scoring_file = parse_options( sys.argv, "scoring_file", "" )
new_torsional_potential= parse_options( sys.argv, "new_torsional_potential", "True" )
cutpoint_open = parse_option_int_list ( sys.argv, 'cutpoint_open' )

if input_pdb =="" :  
    error_exit("USER need to specify -pdb option")
check_path_exist(input_pdb)

if rebuild_res ==0 : 
    error_exit("USER need to specify -rebuild_res option")

if map_file != "" : 
    check_path_exist( map_file )
    map_file = abspath( map_file )

num_pose_kept = 30

native_screen = True 
if native_screen_RMSD > 10.0 :
    native_screen = False
##############location of executable and database:########################
database_folder = rosetta_database() 
rna_swa_test_exe = rosetta_bin("swa_rna_main.linuxgccrelease" )
rna_anal_loop_close_exe = rosetta_bin("swa_rna_analytical_closure.linuxgccrelease" )
#############folder_name##################################################
base_dir=os.getcwd()

main_folder =          os.path.abspath("%s_res_%d" %(basename(input_pdb).replace('.','_'), rebuild_res) ) 
temp_folder=           os.path.abspath("%s/temp_folder/" %(main_folder) )
sampling_folder=       os.path.abspath("%s/sampling/" %(main_folder) ) 
cluster_folder=        os.path.abspath("%s/cluster/" %(main_folder))
output_pdb_folder=     os.path.abspath("%s/output_pdb/" %(main_folder))
precluster_pdb_folder= os.path.abspath("%s/precluster_pdb/" %(main_folder))

if exists(main_folder) : 
    print "warning...main_folder:%s already exist...removing it...! " % main_folder
    remove(main_folder) 
os.mkdir(main_folder)

if exists(temp_folder) : 
    print "warning...temp_folder:%s already exist...removing it...! " % temp_folder
    remove(temp_folder) 
os.mkdir(temp_folder)

if exists(output_pdb_folder) : 
    print "warning...output_pdb_folder:%s already exist...removing it...! " % output_pdb_folder
    remove(output_pdb_folder)
os.mkdir(output_pdb_folder)

if exists(precluster_pdb_folder) : 
    print "warning...precluster_pdb_folder:%s already exist...removing it...! " %(precluster_pdb_folder) 
    remove(precluster_pdb_folder)
os.mkdir(precluster_pdb_folder)

############################################################
input_pdb_temp= temp_folder + '/' + basename(input_pdb)
copy(input_pdb,  input_pdb_temp)
input_pdb=input_pdb_temp

if map_file != "" :
    map_file_temp = temp_folder + '/' + basename(map_file)
    copy(map_file, map_file_temp)
    map_file=map_file_temp

native_pdb = os.path.abspath("%s/native_struct.pdb" % output_pdb_folder)
copy(input_pdb, native_pdb)
#####################Slice out the rebuild region##########################
cutpoint_final = []
res_sliced_all = []
native_pdb_final = native_pdb
rebuild_res_final = rebuild_res
if slice_nearby :
    native_pdb_final = native_pdb.replace('native_struct.pdb', 'native_struct_sliced.pdb')

    rebuilding_res = []
    rebuilding_res.append( rebuild_res )
    if rebuild_res != 1 and (not rebuild_res - 1 in cutpoint_open) :
        rebuilding_res.append( rebuild_res - 1 )
    if rebuild_res != get_total_res(native_pdb) and (not rebuild_res in cutpoint_open) :
        rebuilding_res.append( rebuild_res + 1 )

    res_sliced_all = pdb_slice_with_patching( native_pdb, native_pdb_final, rebuilding_res ) [0]

    rebuild_res_final = res_sliced_all.index(rebuild_res) + 1
    for cutpoint in cutpoint_open :
        if cutpoint in res_sliced_all :
            cutpoint_final.append( res_sliced_all.index(cutpoint) + 1 )
else :
    cutpoint_final = cutpoint_open

total_res = get_total_res(native_pdb_final)
if verbose :
    print "rebuilding_res = %s" % rebuild_res_final
    print "res_sliced = %s" % res_sliced_all
    print "cutpoint_final = %s" % cutpoint_final
    print "total_res = %d " % total_res
#################Check if the rebuilding Rsd is at chain break###########
is_chain_break =False
if rebuild_res_final == 1 or rebuild_res_final == total_res :
    is_chain_break = True
else :
    for cutpoint in cutpoint_final :
        if rebuild_res_final == cutpoint or rebuild_res_final == cutpoint + 1 :
            is_chain_break = True
            break

#Overide ideal_geometry to False when at chain break
if is_chain_break :
    ideal_geometry = False

start_pdb = ""
if ideal_geometry :
    start_pdb = "%s/missing_rebuild_res_native.pdb" % temp_folder 
    pdb_slice(native_pdb_final, start_pdb, "%d-%d %d-%d" % (1, rebuild_res_final-1, rebuild_res_final+1, total_res) )
else :
    start_pdb = native_pdb_final

#####################Create fasta file#####################################
fasta_file=temp_folder + '/fasta'
pdb2fasta(native_pdb_final, fasta_file)
#########################Common Options##################################
common_cmd = "" 

common_cmd += " -database %s " % database_folder

if verbose :
    common_cmd += " -VERBOSE true"
else:
    common_cmd += " -VERBOSE false"

common_cmd += " -fasta %s " % fasta_file

common_cmd += " -input_res "
for n in range(1,total_res+1) : 
    if ideal_geometry and n == rebuild_res_final : 
        continue
    common_cmd += "%d " %n

common_cmd+= " -fixed_res "
for n in range(1,total_res+1) : 
    if n == rebuild_res_final : 
        continue
    common_cmd += "%d " %n

common_cmd += " -jump_point_pairs NOT_ASSERT_IN_FIXED_RES 1-%d " % total_res
common_cmd += " -alignment_res 1-%d " % total_res

common_cmd += " -rmsd_res %d " %(total_res)
common_cmd += " -native " + native_pdb_final
if scoring_file == "" :
    if map_file == "" :
        common_cmd += " -score:weights rna/rna_loop_hires_04092010 "
    else :
        common_cmd += " -score:weights rna/rna_hires_elec_dens "
else :
    common_cmd += " -score:weights %s " % scoring_file


if map_file != "" :
    common_cmd += " -edensity:mapfile %s " % map_file
    common_cmd += " -edensity:mapreso %s " % map_reso
    common_cmd += " -edensity:realign no "

if len(cutpoint_final) != 0 :
    common_cmd += " -cutpoint_open "
    for cutpoint in cutpoint_final :
        common_cmd += '%d ' % cutpoint

if new_torsional_potential :
    common_cmd += " -score:rna_torsion_potential RNA09_based_2012_new "

#########################Sampler Options##################################
if not is_chain_break :
    sampling_cmd = rna_anal_loop_close_exe + ' -algorithm rna_resample_test '
else :
    sampling_cmd = rna_swa_test_exe + ' -algorithm rna_resample_test '

sampling_cmd += " -s %s " % start_pdb
sampling_cmd += " -out:file:silent blah.out " 
sampling_cmd += " -output_virtual true "
sampling_cmd += " -sampler_perform_o2star_pack true "
sampling_cmd += " -sampler_extra_syn_chi_rotamer true "
sampling_cmd += " -sampler_cluster_rmsd %s " % 0.3 
sampling_cmd += " -centroid_screen true "
sampling_cmd += " -minimize_and_score_native_pose %s " % str(include_native).lower()
sampling_cmd += " -finer_sampling_at_chain_closure %s " % str(finer_sampling).lower()
sampling_cmd += " -native_edensity_score_cutoff %s " % native_edensity_cutoff
sampling_cmd += " -sampler_native_rmsd_screen %s " % str(native_screen).lower()
sampling_cmd += " -sampler_native_screen_rmsd_cutoff %s " % native_screen_RMSD
sampling_cmd += " -sampler_num_pose_kept %s " % num_pose_kept 
sampling_cmd += " -PBP_clustering_at_chain_closure true " 
sampling_cmd += " -allow_chain_boundary_jump_partner_right_at_fixed_BP true "
sampling_cmd += " -add_virt_root true "


#############################################################################
if not is_chain_break :
    #################Use Analytical Loop Closure#############################
    if(exists(sampling_folder)): 
        print "warning...sampling_folder:%s already exist...removing it...! " % sampling_folder
        remove(sampling_folder)
    os.mkdir(sampling_folder)

    if is_append :
        print  '\n', "Rebuilding res %s by attaching to res %s" % (rebuild_res_final, rebuild_res_final-1),  '\n'
    else :
        print  '\n', "Rebuilding res %s by attaching to res %s" % (rebuild_res_final, rebuild_res_final+1),  '\n'

    os.chdir( sampling_folder)

    specific_cmd =""
    specific_cmd += " -sample_res %d " % rebuild_res_final
    if is_append :
        specific_cmd += " -cutpoint_closed %d " % rebuild_res_final
    else :
        specific_cmd += " -cutpoint_closed %d " % rebuild_res_final - 1

    command = sampling_cmd + ' ' + specific_cmd + ' ' + common_cmd
    if (verbose): print  '\n', command, '\n'

    subprocess_call( command, 'sampling_1.out', 'sampling_1.err' ) 
     
    os.chdir( base_dir)

###########Sample chainbreak residues with original SWA rebuild###########
else :
    print "Rebuilding residue at chain break point, ignoring chain closure..."
    if(exists(sampling_folder)): 
        print "warning...sampling_folder:%s already exist...removing it...! " % sampling_folder  
        remove(sampling_folder)
    os.mkdir(sampling_folder)

    print  '\n', "Rebuilding res %s" % rebuild_res_final,  '\n'

    os.chdir( sampling_folder)

    specific_cmd=""
    specific_cmd+= " -sample_res %d " % rebuild_res_final

    command = sampling_cmd + ' ' + specific_cmd + ' ' + common_cmd

    if(verbose): print  '\n', command, '\n'

    subprocess_call( command, 'sampling_1.out', 'sampling_1.err' ) 

     
    os.chdir( base_dir)

print '\n',"ALMOST DONE...sorting/clustering/extracting output_pdb....", '\n'
###################Rename the pdb using the clusterer#######################
if(exists(cluster_folder)): 
    print "warning...cluster_folder:%s already exist...removing it...! " % cluster_folder
    remove(cluster_folder)
os.mkdir(cluster_folder)

CONTROL_filename = abspath(output_pdb_folder + "/CONTROL.out")
cluster_filename = abspath(output_pdb_folder + "/cluster.out")
precluster_filename = abspath(precluster_pdb_folder + "/precluster.out")


os.chdir(cluster_folder)

cluster_args = rna_swa_test_exe + " -algorithm rna_cluster "
cluster_args += " -sample_res %d " % rebuild_res_final

if not is_chain_break :
    cluster_args += " -cutpoint_closed %d " % rebuild_res_final

if len(cutpoint_final) != 0 :
    cluster_args+= " -cutpoint_open "
    for cutpoint in cutpoint_final :
        cluster_args += '%d ' % cutpoint

cluster_args += " -rmsd_res %d " % rebuild_res_final
cluster_args += " -add_lead_zero_to_tag true "
cluster_args += " -add_virt_root true "
cluster_args += " -in:file:silent_struct_type  binary_rna"
cluster_args += " -in:file:silent %s/blah.out " % sampling_folder
cluster_args += " -PBP_clustering_at_chain_closure true " 
cluster_args += " -allow_chain_boundary_jump_partner_right_at_fixed_BP true "

no_clustering  = " -suite_cluster_radius 0.0 " 
no_clustering += " -loop_cluster_radius 0.0 "

if(verbose):  ##This is just for control purposes...
    command = (cluster_args + ' ' + common_cmd + no_clustering + 
      " -recreate_silent_struct false  -out:file:silent %s" % CONTROL_filename)

    if exists(sampling_folder + '/blah.out'):
        print '\n', command ,'\n'
        subprocess_call( command, 'CONTROL.out', 'CONTROL.err' )

    ###This one is with alignment to native_pose...however wary that SCORE output is not correct...
    command = (cluster_args + ' ' + common_cmd +  no_clustering + 
      " -recreate_silent_struct true -out:file:silent %s" % precluster_filename)

    if exists(sampling_folder + '/blah.out'):
        if(verbose): print '\n', command ,'\n'
        subprocess_call( command, 'precluster.out', 'precluster.err' )

with_clustering  = ""
with_clustering += " -suite_cluster_radius %s " % cluster_RMSD
with_clustering += " -loop_cluster_radius 999.99 "
with_clustering += " -clusterer_num_pose_kept 50 "

command = (cluster_args + ' ' + common_cmd +  with_clustering + 
  " -recreate_silent_struct true -out:file:silent %s" % cluster_filename )

if exists(sampling_folder + '/blah.out'):
    if(verbose): print '\n', command ,'\n'
    subprocess_call( command, 'cluster.out', 'cluster.err' )

os.chdir( base_dir)

if verbose : 
    if exists(CONTROL_filename) :
        score_line = subprocess_out('head -n 2 %s ' % CONTROL_filename ) [1]
        subprocess_call('echo "%s" > %s/output_pdb_CONTROL.txt ' % (score_line, main_folder))
        subprocess_call("grep SCORE %s | sort -nk2 >> %s/output_pdb_CONTROL.txt " % (CONTROL_filename, main_folder))

    if exists(precluster_filename) :
        score_line = subprocess_out('head -n 2 %s ' % precluster_filename ) [1]
        subprocess_call('echo "%s" > %s/output_pdb_precluster.txt '  %(score_line, main_folder))
        subprocess_call("grep SCORE %s | sort -nk2 >> %s/output_pdb_precluster.txt " % (precluster_filename, main_folder))

if exists(cluster_filename) :
    score_line = subprocess_out('head -n 2 %s ' % cluster_filename ) [1]
    subprocess_call('echo "%s" > %s/output_pdb.txt ' % (score_line, main_folder) )
    subprocess_call("grep SCORE %s | sort -nk2 >> %s/output_pdb.txt " % (cluster_filename, main_folder))
else :
    output = open("%s/output_pdb.txt" % main_folder, 'w')
    output.write("No silent file is being output during rebuilding")
#####################Extract the pdb from the silent_file###################
if verbose and exists(precluster_filename) :
    extract_pdb(precluster_filename, precluster_pdb_folder)
if exists(cluster_filename) :
    extract_pdb(cluster_filename, output_pdb_folder)
#####################Merge the sliced region back to starting pdb##################
if slice_nearby :
    os.chdir( output_pdb_folder )
    pdb_file_list = glob("*.pdb")
    for pdb_file in pdb_file_list :
        sliced2orig_merge_back( native_pdb, pdb_file, pdb_file.replace('.pdb', '_merge.pdb'), res_sliced_all )
    os.chdir(base_dir)
########################################################################

if not verbose :
    if exists(cluster_filename):
        remove(cluster_filename)
    remove(temp_folder)
    remove(sampling_folder)
    remove(cluster_folder)
    remove(precluster_pdb_folder)

total_time=time.time()-start_time

print '\n', "DONE!...Total time taken= %f seconds" % total_time , '\n'     
print '###################################'	


