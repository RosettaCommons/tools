#!/usr/bin/env python
import os.path
import imp

try :
    import erraser_util
except :
    file_path = os.path.split( os.path.abspath(__file__) ) [0]
    imp.load_source('erraser_util', file_path + '/erraser_util.py')

from erraser_util import *

#####################################################
print '###################################'	
print 'Starting SWA_rebuild_erraser.py...'	
start_time=time.time()


#####Input options#####################################
input_pdb = parse_options(sys.argv, 'pdb', '')
out_pdb = parse_options(sys.argv, 'out_pdb', basename(input_pdb).replace('.pdb', '_erraser.pdb') )
n_iterate = parse_options(sys.argv, 'n_iterate', 1)
map_file = parse_options(sys.argv, 'map', '')
map_reso = parse_options(sys.argv, 'map_reso', '1.5')
native_screen_RMSD= parse_options(sys.argv, "native_screen_RMSD", 2.0)
finer_sampling = parse_options(sys.argv, 'finer_sampling', 'False')
new_torsional_potential = parse_options( sys.argv, "new_torsional_potential", "True" )
use_existing_temp_folder = parse_options( sys.argv, "use_existing_temp_folder", "True" )
kept_temp_folder = parse_options( sys.argv, "kept_temp_folder", "False" )
rebuild_all = parse_options( sys.argv, "rebuild_all", "False" )
rebuild_rmsd = parse_options( sys.argv, "rebuild_rmsd", "True" )
extra_res = parse_option_chain_res_list ( sys.argv, 'rebuild_extra_res' )
fixed_res = parse_option_chain_res_list ( sys.argv, 'fixed_res' )
cutpoint_open = parse_option_chain_res_list ( sys.argv, 'cutpoint_open' )

if input_pdb =="" : 
    error_exit("Error: USER need to specify -pdb option")
if map_file =="" : 
    error_exit("Error: USER need to specify -map option")
check_path_exist(input_pdb)
check_path_exist(map_file)

if exists( out_pdb ) : 
    print 'Output pdb %s already exists... Remove it.' % out_pdb
    remove(out_pdb)
#######File Paths#######################################################
python_file_path = os.path.split( os.path.abspath(__file__) ) [0]
minimize_python = "%s/full_struct_slice_and_minimize.py" % python_file_path
seq_rebuild_python = "%s/seq_rebuild.py" % python_file_path
check_path_exist(minimize_python)
check_path_exist(seq_rebuild_python)
input_pdb = abspath(input_pdb)
map_file = abspath(map_file)
out_pdb = abspath(out_pdb)
#####Temporary data output folder#######################################
base_dir = os.getcwd()
temp_dir = '%s/%s_erraser_temp/' % (base_dir, basename(input_pdb).replace('.pdb', ''))
if exists(temp_dir) :
    if use_existing_temp_folder :
        print 'Temporay directory %s exists... Use the data stored in the existing folder.' % temp_dir
        print 'Because -use_existing_temp_folder is set to True.'
    else :
        print 'Temporay directory %s exists... Remove it and create a new folder.' % temp_dir
        print 'Because -use_existing_temp_folder is set to False.'
        remove(temp_dir)
        os.mkdir(temp_dir)
else :
    print 'Create temporay directory %s...' % temp_dir
    os.mkdir(temp_dir) 
########################################################################
os.chdir(temp_dir)
regularized_input_pdb = basename(input_pdb).replace('.pdb', '_regularized.pdb')
regularize_pdb(input_pdb, regularized_input_pdb)
[res_conversion_list, fixed_res_final, cutpoint_final, CRYST1_line] = pdb2rosetta(regularized_input_pdb, 'start.pdb')
if not exists('minimize_0.pdb') :
    rna_rosetta_ready_set('start.pdb', 'minimize_0.pdb')
else :
    print 'minimize_0.pdb already exists... Skip the ready_set step.'

for res in fixed_res :
    if res in res_conversion_list :
        fixed_res_final.append( res_conversion_list.index(res) + 1 )

for res in cutpoint_open :
    if res in res_conversion_list :
        cutpoint_final.append( res_conversion_list.index(res) + 1 )

extra_res_final = []
for res in extra_res :
    if res in res_conversion_list :
        extra_res_final.append( res_conversion_list.index(res) + 1 )

fixed_res_final.sort()
cutpoint_final.sort()
extra_res_final.sort()
print "fixed_res_final in Rosetta pdb file = %s" % fixed_res_final
print "cutpoint_final in Rosetta pdb file = %s" % cutpoint_final
print "extra_res_final in Rosetta pdb file = %s" % extra_res_final

for res in extra_res_final :
    if res in fixed_res_final :
        res_name = res_conversion_list[res - 1]
        error_message  = "Confliction: rebuild_extra_res %s is either covalently bonded to a modified base" % res_name
        error_message += " or in user-input -fixed_res!!!" 
        error_exit(error_message)
#############################################################################
minimize_command_common  = minimize_python
minimize_command_common += " -map %s" % map_file
minimize_command_common += " -map_reso %s" % map_reso
minimize_command_common += ' -new_torsional_potential %s' % new_torsional_potential
if len(fixed_res_final) != 0 :
    minimize_command_common += " -fixed_res"
    for res in fixed_res_final :
        minimize_command_common += " %s" % res
        
minimize_command  = minimize_command_common
minimize_command += " -pdb minimize_0.pdb"
minimize_command += " -out_pdb minimize_1.pdb"

if not exists('minimize_1.pdb') :
    print 'Starting whole-structure minimization 1...' 
    subprocess_call(minimize_command, 'minimize_1.out', 'minimize_1.err')
else :
    print 'minimize_1.pdb already exists... Skip the minimization step.'
print 'Minimization 1 completed sucessfully.'
#############################################################################

for step in range(1, n_iterate + 1) :

    #Find errorenous nucleotides for rebuilding using phenix.rna_validate
    rosetta2std_pdb( 'minimize_%s.pdb' % step, 'standardized.pdb', CRYST1_line )
    rebuild_res_error = find_error_res("standardized.pdb")
    #Also include extra_res in the error category
    for res in extra_res_final :
        if not res in rebuild_res_error :
            rebuild_res_error.append(res)

    #Find residues with large RMSD before and after minimization for rebuilding
    #Excludes residues already in rebuild_res_error
    rebuild_res_rmsd = []
    if rebuild_all :
        ###Overide the RMSD selection and rebuild all residues if rebuild_all == True###
        total_res = get_total_res('minimize_%s.pdb' % step)
        rebuild_res_rmsd = range(1, total_res + 1)
    elif rebuild_rmsd :
        res_rmsd_list = res_wise_rmsd('minimize_%s.pdb' % (step - 1),'minimize_%s.pdb' % step)
        res_rmsd_list_temp = []
        for res_list in res_rmsd_list :
            if not res_list[0] in rebuild_res_error :
                res_rmsd_list_temp.append(res_list)
        res_rmsd_list = sorted(res_rmsd_list_temp, key = lambda x : x[1], reverse=True)
        rmsd_cutoff = float(map_reso) * 0.05
        percentage_cutoff = 0.2

        for i in range(0, int(len(res_rmsd_list) * percentage_cutoff)) :
            if res_rmsd_list[i][1] > rmsd_cutoff :
                rebuild_res_rmsd.append(res_rmsd_list[i][0])

    #Residues in fixed_res are not rebuilt so needed to be removed from the lists
    for res in fixed_res_final :
        if res in rebuild_res_error :
            rebuild_res_error.remove(res)
        if res in rebuild_res_rmsd :
            rebuild_res_rmsd.remove(res)
    #Remove residues in rebuild_res_rmsd that are overlaped with rebuild_res_error
    for res in rebuild_res_error :
        if res in rebuild_res_rmsd :
            rebuild_res_rmsd.remove(res)
    rebuild_res_error.sort()
    rebuild_res_rmsd.sort()

    if len(rebuild_res_error) == 0 and len(rebuild_res_rmsd) == 0 :
        print "No residue need to be rebuilt!"
        copy('minimize_%s.pdb' % step, 'FINAL.pdb')
        rosetta2phenix_merge_back(regularized_input_pdb, 'FINAL.pdb', out_pdb)
        if not kept_temp_folder :
            os.chdir(base_dir)
            remove(temp_dir)
        total_time=time.time()-start_time
        print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
        print '###################################'	
        sys.exit(0)
#############################################################################
    seq_rebuild_command_common = seq_rebuild_python
    seq_rebuild_command_common += " -map %s" % map_file
    seq_rebuild_command_common += " -map_reso %s" % map_reso
    seq_rebuild_command_common += " -finer_sampling %s" % finer_sampling
    seq_rebuild_command_common += ' -new_torsional_potential %s' % new_torsional_potential
    seq_rebuild_command_common += ' -native_screen_RMSD %s' % native_screen_RMSD

    if len(cutpoint_final) != 0 :
        seq_rebuild_command_common += " -cutpoint_open"
        for res in cutpoint_final :
            seq_rebuild_command_common += " %s" % res

    seq_rebuild_command1  = seq_rebuild_command_common
    seq_rebuild_command1 += " -pdb minimize_%s.pdb" % step
    seq_rebuild_command1 += " -out_pdb rebuild_outlier_%s.pdb" % step
    seq_rebuild_command1 += " -rebuild_res"
    for res in rebuild_res_error :
        seq_rebuild_command1 += " %s" % res

    seq_rebuild_command2  = seq_rebuild_command_common
    seq_rebuild_command2 += " -pdb rebuild_outlier_%s.pdb" % step
    seq_rebuild_command2 += " -out_pdb rebuild_%s.pdb" % step
    seq_rebuild_command2 += " -native_edensity_cutoff 0.97"
    seq_rebuild_command2 += " -rebuild_res"
    for res in rebuild_res_rmsd :
        seq_rebuild_command2 += " %s" % res

    if not exists('rebuild_%s.pdb' % step) :
        if not exists('rebuild_outlier_%s.pdb' % step) :

            if len(rebuild_res_error) != 0 :
                subprocess_call(seq_rebuild_command1, 'rebuild_outlier_%s.out' % step, 'rebuild_outlier_%s.err' % step)
            else :
                print 'No errorenous residues... Skip the error residues rebuilding step %s.' % step
                copy('minimize_%s.pdb' % step, 'rebuild_outlier_%s.pdb' % step)

        else :
            print 'rebuild_outlier_%s.pdb already exists... Skip outlier rebuilding step.' % step

        if len(rebuild_res_rmsd) != 0 :
            subprocess_call(seq_rebuild_command2, 'rebuild_rmsd_%s.out' % step, 'rebuild_rmsd_%s.err' % step)
        else :
            if rebuild_rmsd :
                print 'No high-RMSD residues... Skip the high-RMSD residues rebuilding step %s.' % step
            else :
                print 'rebuild_rmsd=False... Skip the high-RMSD residues rebuilding step %s.' % step
            copy('rebuild_outlier_%s.pdb' % step, 'rebuild_%s.pdb' % step)

    else :
        print 'rebuild_%s.pdb already exists... Skip rebuilding step.' % step
    print 'Rebuilding %s completed sucessfully.' % step
##############################################################################
    minimize_command  = minimize_command_common
    minimize_command += " -pdb rebuild_%s.pdb" % step
    minimize_command += " -out_pdb minimize_%s.pdb" % (step +1)

    if not exists( 'minimize_%s.pdb' % (step + 1) ) :
        subprocess_call( minimize_command, 'minimize_%s.out' % (step + 1), 'minimize_%s.err' % (step + 1) )
    else :
        print 'minimize_%s.pdb already exists... Skip the minimization step.' % (step + 1)
    print 'Minimization %s completed sucessfully.' % (step + 1)
##############################################################################

copy( 'minimize_%s.pdb' % (step + 1), 'FINAL.pdb' )
rosetta2phenix_merge_back(regularized_input_pdb, 'FINAL.pdb', out_pdb)
if not kept_temp_folder :
    os.chdir(base_dir)
    remove(temp_dir)

total_time=time.time()-start_time
print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
print '###################################'	

