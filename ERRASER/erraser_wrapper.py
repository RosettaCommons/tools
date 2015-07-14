import os.path
import imp
from copy import deepcopy

file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('erraser_util', file_path + '/erraser_util.py')
imp.load_source('erraser_option', file_path + '/erraser_option.py')
from erraser_util import *


##### erraser start####################################################
def erraser( option ) :

    stdout = sys.stdout
    stderr = sys.stderr
    if option.log_out != "" :
        sys.stdout = open(option.log_out, 'w')
    if option.log_err != "" :
        sys.stderr = open(option.log_err, 'w')

    #####################################################
    print '###################################'
    print 'Starting erraser...'
    start_time=time.time()
    #Faster rebuilding mode
    option.num_pose_kept_cluster = 1
    option.finalize()

    #####Temporary data output folder##################################
    base_dir = os.getcwd()
    temp_dir = '%s/%s_erraser_temp/' % (base_dir, basename(option.input_pdb).replace('.pdb', ''))
    if exists(temp_dir) :
        if option.use_existing_temp_folder :
            print 'Temporary directory %s exists... Use the data stored in the existing folder.' % temp_dir
            print 'Because -use_existing_temp_folder is set to True.'
        else :
            print 'Temporary directory %s exists... Remove it and create a new folder.' % temp_dir
            print 'Because -use_existing_temp_folder is set to False.'
            remove(temp_dir)
            os.mkdir(temp_dir)
    else :
        print 'Create temporary directory %s...' % temp_dir
        os.mkdir(temp_dir)
    ########################################################################
    os.chdir(temp_dir)
    regularized_input_pdb = basename(option.input_pdb).replace('.pdb', '_regularized.pdb')
    regularize_pdb(option.input_pdb, regularized_input_pdb)
    [res_conversion_list, fixed_res_final, cutpoint_final, CRYST1_line] = pdb2rosetta(regularized_input_pdb, 'start.pdb', using_protein = option.rna_prot_erraser)
    if not exists('minimize_0.pdb') :
        rna_rosetta_ready_set('start.pdb', 'minimize_0.pdb', option.rosetta_bin, option.rosetta_database, option.rna_prot_erraser)
    else :
        print 'minimize_0.pdb already exists... Skip the ready_set step.'

    for res in option.fixed_res :
        fixed_res_final.append( res_num_convert(res_conversion_list, res) )

    for res in option.cutpoint :
        cutpoint_final.append( res_num_convert(res_conversion_list, res) )

    extra_res_final = []
    for res in option.extra_res :
        extra_res_final.append( res_num_convert(res_conversion_list, res) )

    fixed_res_final.sort()
    cutpoint_final.sort()
    extra_res_final.sort()
    print "fixed_res_final in Rosetta pdb file = %s" % fixed_res_final
    print "cutpoint_final in Rosetta pdb file = %s" % cutpoint_final
    print "extra_res_final in Rosetta pdb file = %s" % extra_res_final
    option.fixed_res_rs = fixed_res_final
    option.cutpoint_rs = cutpoint_final

    for res in extra_res_final :
        if res in fixed_res_final :
            res_name = res_conversion_list[res - 1]
            error_message  = "Confliction: rebuild_extra_res %s is either covalently bonded to a modified base" % res_name
            error_message += " or in user-input -fixed_res!!!"
            error_exit(error_message)
    ###################################################################
    minimize_option = deepcopy(option)
    rebuild_option = deepcopy(option)
    rebuild_option.ideal_geometry = False
    if not exists('minimize_1.pdb') :
        print 'Starting whole-structure minimization 1...'
        minimize_option.input_pdb = "minimize_0.pdb"
        minimize_option.out_pdb = "minimize_1.pdb"
        minimize_option.log_out = "minimize_1.out"
        full_struct_slice_and_minimize( minimize_option )
    else :
        print 'minimize_1.pdb already exists... Skip the minimization step.'
    print 'Minimization 1 completed sucessfully.'
    ###################################################################

    for step in range(1, option.n_iterate + 1) :

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
        if option.rebuild_all :
            ###Overide the RMSD selection and rebuild all residues if rebuild_all == True###
            total_res = get_total_res('minimize_%s.pdb' % step)
            rebuild_res_rmsd = range(1, total_res + 1)
        elif option.rebuild_rmsd :
            res_rmsd_list = res_wise_rmsd('minimize_%s.pdb' % (step - 1),'minimize_%s.pdb' % step)
            res_rmsd_list_temp = []
            for res_list in res_rmsd_list :
                if not res_list[0] in rebuild_res_error :
                    res_rmsd_list_temp.append(res_list)
            res_rmsd_list = sorted(res_rmsd_list_temp, key = lambda x : x[1], reverse=True)
            rmsd_cutoff = float(option.map_reso) * 0.05
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
            rosetta2phenix_merge_back(regularized_input_pdb, 'FINAL.pdb', option.out_pdb)
            os.chdir(base_dir)
            if not option.kept_temp_folder :
                remove(temp_dir)
            total_time=time.time()-start_time
            print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
            print '###################################'
            sys.exit(0)
        ###################################################################
        if not exists('rebuild_%s.pdb' % step) :
            if not exists('rebuild_outlier_%s.pdb' % step) :

                if len(rebuild_res_error) != 0 :
                    rebuild_option.input_pdb = "minimize_%s.pdb" % step
                    rebuild_option.out_pdb = "rebuild_outlier_%s.pdb" % step
                    rebuild_option.rebuild_res_list = rebuild_res_error
                    rebuild_option.log_out = 'rebuild_outlier_%s.out' % step
                    seq_rebuild( rebuild_option )
                else :
                    print 'No errorenous residues... Skip the error residues rebuilding step %s.' % step
                    copy('minimize_%s.pdb' % step, 'rebuild_outlier_%s.pdb' % step)

            else :
                print 'rebuild_outlier_%s.pdb already exists... Skip outlier rebuilding step.' % step

            if len(rebuild_res_rmsd) != 0 :
                rebuild_option.input_pdb = "rebuild_outlier_%s.pdb" % step
                rebuild_option.out_pdb = "rebuild_%s.pdb" % step
                rebuild_option.rebuild_res_list = rebuild_res_rmsd
                rebuild_option.native_edensity_cutoff = 0.97
                rebuild_option.log_out = 'rebuild_rmsd_%s.out' % step
                seq_rebuild( rebuild_option )
            else :
                if option.rebuild_rmsd :
                    print 'No high-RMSD residues... Skip the high-RMSD residues rebuilding step %s.' % step
                else :
                    print 'rebuild_rmsd=False... Skip the high-RMSD residues rebuilding step %s.' % step
                copy('rebuild_outlier_%s.pdb' % step, 'rebuild_%s.pdb' % step)

        else :
            print 'rebuild_%s.pdb already exists... Skip rebuilding step.' % step
        print 'Rebuilding %s completed sucessfully.' % step
        ###################################################################
        if not exists( 'minimize_%s.pdb' % (step + 1) ) :
            minimize_option.input_pdb = "rebuild_%s.pdb" % step
            minimize_option.out_pdb = "minimize_%s.pdb" % (step +1)
            minimize_option.log_out = "minimize_%s.out" % (step +1)
            full_struct_slice_and_minimize( minimize_option )
        else :
            print 'minimize_%s.pdb already exists... Skip the minimization step.' % (step + 1)
        print 'Minimization %s completed sucessfully.' % (step + 1)
    ###################################################################

    copy( 'minimize_%s.pdb' % (step + 1), 'FINAL.pdb' )
    rosetta2phenix_merge_back(regularized_input_pdb, 'FINAL.pdb', 'FINAL_merge.pdb')
    regularize_OP1_OP2('FINAL_merge.pdb', option.out_pdb)
    os.chdir(base_dir)
    if not option.kept_temp_folder :
        remove(temp_dir)

    total_time=time.time()-start_time
    print '\n', "DONE!...Total time taken= %f seconds" %(total_time) , '\n'
    print '###################################'
    if sys.stdout != sys.__stdout__:
        sys.stdout.close()
    if sys.stderr != sys.__stderr__:
        sys.stderr.close()
    sys.stdout = stdout
    sys.stderr = stderr
##### erraser end #####################################################



##### erraser_single_res start ########################################
def erraser_single_res( option ) :

    stdout = sys.stdout
    stderr = sys.stderr
    if option.log_out != "" :
        sys.stdout = open(option.log_out, 'w')
    if option.log_err != "" :
        option.sys.stderr = open(option.log_err, 'w')
    #Search all Chi conformer without constraint
    option.search_syn_pyrimidine_only_when_native_syn = False
    option.constrain_chi = False

    print '###################################'
    print 'Starting erraser_single_res...'
    start_time=time.time()

    option.finalize()
    if option.out_prefix == "" :
        option.out_prefix = basename(option.input_pdb).replace('.pdb', '_rebuild')
    option.out_prefix = abspath(option.out_prefix)
    #####Set temp folder#######################
    base_dir = os.getcwd()
    temp_dir = '%s/%s/' % (base_dir, basename(option.input_pdb).replace('.pdb', '_single_res_rebuild_temp') )
    if exists(temp_dir) :
        print 'Temporary directory %s exists... Remove it and create a new folder.' % temp_dir
        remove(temp_dir)
        os.mkdir(temp_dir)
    else :
        print 'Create temporary directory %s...' % temp_dir
        os.mkdir(temp_dir)

    print '###################################'
    #####################################################
    os.chdir(temp_dir)

    regularized_input_pdb = basename(option.input_pdb).replace('.pdb', '_regularized.pdb')
    regularize_pdb(option.input_pdb, regularized_input_pdb)
    [res_conversion_list, fixed_res_final, cutpoint_final, CRYST1_line] = pdb2rosetta(regularized_input_pdb, 'start.pdb', using_protein = option.rna_prot_erraser)
    cutpoint_final.sort()
    fixed_res_final.sort()
    option.fixed_res_rs = fixed_res_final
    option.cutpoint_rs = cutpoint_final
    option.rebuild_res = res_num_convert(res_conversion_list, option.rebuild_res_pdb)
    rna_rosetta_ready_set('start.pdb', 'temp.pdb', option.rosetta_bin, option.rosetta_database)

    print 'Starting to rebuild residue %s' % option.rebuild_res_pdb

    SWA_option = deepcopy(option)
    SWA_option.input_pdb = 'temp.pdb'
    SWA_option.log_out = 'rebuild_%s.out' % option.rebuild_res_pdb
    SWA_rebuild_erraser( SWA_option )
    rebuilt_pdb_merge0 = './temp_pdb_res_%d/output_pdb/S_000000_merge.pdb' % option.rebuild_res
    rebuilt_pdb_orig0 = './temp_pdb_res_%d/output_pdb/S_000000.pdb' % option.rebuild_res
    rebuilt_pdb_final = ''

    if exists(rebuilt_pdb_merge0) :
        rebuilt_pdb_final = rebuilt_pdb_merge0
    elif exists(rebuilt_pdb_orig0) :
        rebuilt_pdb_final = rebuilt_pdb_orig0
    else :
        print "No alternative conformation is found...."

    if rebuilt_pdb_final != '' :
        print "Residue %s is sucessfully rebuilt!" % option.rebuild_res_pdb
        res_sliced = []
        rebuilding_res_rs = 0
        for line in open(SWA_option.log_out) :
            if "res_sliced" in line :
                for elem in line.split('=') [-1].split(',') :
                    elem = elem.replace('[', '').replace(']', '')
                    res_sliced.append( int(elem) )
            if "rebuilding_res" in line :
                rebuilding_res_rs = int( line.split() [-1] )

        fixed_res = []
        total_res = get_total_res(rebuilt_pdb_orig0)
        if rebuilding_res_rs == 1 :
            fixed_res = range(2,total_res+1)
        elif rebuilding_res_rs == total_res :
            fixed_res = range(1, total_res)
        else :
            fixed_res = range(1, rebuilding_res_rs) + range(rebuilding_res_rs+1, total_res+1)
        option.fixed_res_rs = fixed_res

        final_pdb_list = []
        native_pdb = './temp_pdb_res_%d/output_pdb/native_struct.pdb' % option.rebuild_res
        native_pdb_sliced = './temp_pdb_res_%d/output_pdb/native_struct_sliced.pdb' % option.rebuild_res

        minimize_option = deepcopy(option)

        #Scores for other rebuilt models
        for i in xrange(option.num_pose_kept_cluster) :
            rebuilt_pdb = './temp_pdb_res_%d/output_pdb/S_%06d.pdb' % (option.rebuild_res, i)
            if exists(rebuilt_pdb) :
                if option.skip_single_res_minimize:
                    if len(res_sliced) != 0 :
                        merged_pdb = rebuilt_pdb.replace('.pdb', '_merge.pdb')
                        sliced2orig_merge_back(native_pdb, rebuilt_pdb, merged_pdb, res_sliced )
                    final_pdb_list.append([i, merged_pdb])
                else:
                    minimize_option.input_pdb = rebuilt_pdb
                    minimize_option.out_pdb = rebuilt_pdb.replace('.pdb', 'min.pdb')
                    minimize_option.log_out = './temp_pdb_res_%d/output_pdb/minimize_%d.out' % (option.rebuild_res, i)
                    erraser_minimize( minimize_option )
                    score = 0.0
                    min_out_lines = open(minimize_option.log_out).readlines()
                    for j in xrange( len(min_out_lines) - 1, -1, -1) :
                        if "current_score =" in min_out_lines[j] or "Total weighted score:" in min_out_lines[j]:
                            score = float(min_out_lines[j].split()[-1])
                            break
                    if len(res_sliced) != 0 :
                        sliced2orig_merge_back( native_pdb, minimize_option.out_pdb, rebuilt_pdb.replace('.pdb', '_merge.pdb'), res_sliced )
                        final_pdb_list.append([score, rebuilt_pdb.replace('.pdb', '_merge.pdb')])
                    else :
                        final_pdb_list.append([score, minimize_option.out_pdb])
            else :
                break
        final_pdb_list = sorted(final_pdb_list)

        #Get (minimized) native score
        if option.skip_single_res_minimize:
            native_merge_pdb = native_pdb
            native_score = 0.0
        else:
            native_merge_pdb = ''
            if exists(native_pdb_sliced) :
                minimize_option.input_pdb = native_pdb_sliced
            else :
                minimize_option.input_pdb = native_pdb
            minimize_option.out_pdb = native_pdb.replace('.pdb', 'min.pdb')
            minimize_option.log_out = './temp_pdb_res_%d/output_pdb/minimize_native.out' % option.rebuild_res
            erraser_minimize( minimize_option )
            native_score = 0.0
            min_out_lines = open(minimize_option.log_out).readlines()
            for j in xrange( len(min_out_lines) - 1, -1, -1) :
                if "current_score =" in min_out_lines[j] or "Total weighted score:" in min_out_lines[j]:
                    native_score = float(min_out_lines[j].split()[-1])
                    break
            if len(res_sliced) != 0 :
                sliced2orig_merge_back( native_pdb, minimize_option.out_pdb, native_pdb.replace('.pdb', '_merge.pdb'), res_sliced )
                native_merge_pdb = native_pdb.replace('.pdb', '_merge.pdb')
            else :
                native_merge_pdb = minimize_option.out_pdb

        #Output scores
        out_score = open("../scores.out" ,'w')
        for i in xrange(0, len(final_pdb_list) ) :
            out_score.write("model_%d %.3f\n" % (i, final_pdb_list[i][0]) )
            rosetta2phenix_merge_back(regularized_input_pdb, final_pdb_list[i][1], "FINAL_merge.pdb")
            regularize_OP1_OP2('FINAL_merge.pdb',  "%s_%d.pdb" % (option.out_prefix, i))
        out_score.write("start_min %.3f\n" % native_score )
        rosetta2phenix_merge_back(regularized_input_pdb, native_merge_pdb, "FINAL_merge.pdb")
        regularize_OP1_OP2('FINAL_merge.pdb',  "%s_start_min.pdb" % option.out_prefix)
        out_score.close()

    if not option.verbose :
        remove('temp_pdb_res_%d' % option.rebuild_res)

    print '###################################'

    os.chdir(base_dir)

    if not option.kept_temp_folder :
        remove(temp_dir)

    total_time=time.time()-start_time
    print '\n', "DONE!...Total time taken= %f seconds" % total_time
    print '###################################'
    if sys.stdout != sys.__stdout__:
        sys.stdout.close()
    if sys.stderr != sys.__stderr__:
        sys.stderr.close()
    sys.stdout = stdout
    sys.stderr = stderr
##### erraser_single_res end   ########################################




##### erraser_minimize start ##########################################
def erraser_minimize( option ) :

    stdout = sys.stdout
    stderr = sys.stderr
    if option.log_out != "" :
        sys.stdout = open(option.log_out, 'w')
    if option.log_err != "" :
        sys.stderr = open(option.log_err, 'w')

    print '###################################'
    print 'Starting erraser_minimize...'
    start_time=time.time()
    option.finalize()

    #######Folders and files paths###########################
    database_folder = rosetta_database_path(option.rosetta_database)
    rna_minimize_exe = rosetta_bin_path("erraser_minimizer", option.rosetta_bin)
    temp_rs = option.input_pdb.replace('.pdb', '_temp_rs.pdb')
    temp_rs_min = option.input_pdb.replace('.pdb', '_temp_rs_min.pdb')

    if exists(temp_rs) :
        print "Temporary file %s exists... Remove it..." % temp_rs
        remove(temp_rs)
    if exists(temp_rs_min) :
        print "Temporary file %s exists... Remove it..." % temp_rs_min
        remove(temp_rs_min)
    ####slicing into smaller pdbs if given in the option####
    fixed_res_final = []
    res_sliced_all = []
    if len(option.res_slice) != 0 :
        [res_sliced_all, patched_res] = pdb_slice_with_patching( option.input_pdb, temp_rs, option.res_slice )
        fixed_res_final += patched_res
        for res in option.fixed_res_rs :
            if res in res_sliced_all :
                fixed_res_final.append( res_sliced_all.index(res) + 1 )
    else :
        copy(option.input_pdb, temp_rs)
        fixed_res_final = option.fixed_res_rs
    option.fixed_res_rs = fixed_res_final
    ####submit rosetta cmdline##############
    command = rna_minimize_exe
    command += " -database %s " % database_folder
    command += " -native %s " % temp_rs
    command += " -out_pdb %s " % temp_rs_min
    command += " -score:weights %s " % option.scoring_file

    if option.fcc2012_new_torsional_potential :
        command += " -score:rna_torsion_potential FCC2012_RNA11_based_new "
    elif option.new_torsional_potential :
        command += " -score:rna_torsion_potential RNA11_based_new "

    command += " -rna::corrected_geo %s " % str(option.corrected_geo).lower()
    command += " -rna::rna_prot_erraser %s " % str(option.rna_prot_erraser).lower()
    command += " -vary_geometry %s " % str(option.vary_geometry).lower()
    command += " -constrain_P %s " % str(option.constrain_phosphate).lower()
    command += " -attempt_pyrimidine_flip %s " % str(option.attempt_pyrimidine_flip).lower()
    command += " -skip_minimize %s " % str(option.skip_minimize).lower()

    #Rescue the minimization default before r53221
    command += " -scale_d 100 "
    command += " -scale_theta 10 "

    #Rescue 2012 defaults
    if option.enlarge_H_lj is True:
        command += " -enlarge_H_lj %s " % str(option.enlarge_H_lj).lower()
    if option.enlarge_lj_hbond_dis is True:
        command += " -lj_hbond_hdis 1.95 "
        command += " -lj_hbond_OH_donor_dis 3.0 "
    if option.restore_pre_talaris_2013_behavior is True:
        command += " -restore_pre_talaris_2013_behavior %s " % str(option.restore_pre_talaris_2013_behavior).lower()
        command += " -score::smooth_fa_elec %s " % str(option.smooth_fa_elec).lower()
    if option.use_2prime_OH_potential is False:
        command += " -use_2prime_OH_potential %s " % str(option.use_2prime_OH_potential).lower()

    if len(option.fixed_res_rs) != 0 :
        command += ' -fixed_res '
        for i in option.fixed_res_rs :
            command += '%d ' % i

    if option.map_file != '' :
        command += " -edensity:mapfile %s " % option.map_file
        command += " -edensity:mapreso %s " % option.map_reso
        command += " -edensity:realign no "
    print "cmdline: %s" % command
    print "#######Submit the Rosetta Command Line###############"
    subprocess_call(command, sys.stdout, sys.stderr)
    print "#####################################################"
    ####Merge final result back to pdb####
    if len(res_sliced_all) != 0 :
        sliced2orig_merge_back( option.input_pdb, temp_rs_min, option.out_pdb, res_sliced_all )
        remove( temp_rs_min )
    else :
        move(temp_rs_min, option.out_pdb)
    remove( temp_rs )
    #########################################
    total_time=time.time()-start_time
    print '\n', "DONE!...Total time taken= %f seconds" % total_time
    print '###################################'
    if sys.stdout != sys.__stdout__:
        sys.stdout.close()
    if sys.stderr != sys.__stderr__:
        sys.stderr.close()
    sys.stdout = stdout
    sys.stderr = stderr
##### erraser_minimize end   ##########################################




##### full_struct_slice_and_minimize start ############################
def full_struct_slice_and_minimize( option ) :

    stdout = sys.stdout
    stderr = sys.stderr
    if option.log_out != "" :
        sys.stdout = open(option.log_out, 'w')
    if option.log_err != "" :
        sys.stderr = open(option.log_err, 'w')

    print '###################################'
    print 'Starting full_struct_slice_and_minimize...'
    start_time=time.time()
    option.finalize()

    #####Set temp folder#######################
    base_dir = os.getcwd()
    temp_dir = '%s/%s' % (base_dir, basename(option.input_pdb).replace('.pdb', '_full_minimize_temp') )
    if exists(temp_dir) :
        print 'Temporary directory %s exists... Remove it and create a new folder.' % temp_dir
        remove(temp_dir)
        os.mkdir(temp_dir)
    else :
        print 'Create temporary directory %s...' % temp_dir
        os.mkdir(temp_dir)
    ###########################################
    os.chdir(temp_dir)
    copy( option.input_pdb, 'before_min.pdb' )
    total_res = get_total_res(option.input_pdb)
    n_chunk = int(total_res / 100.0 + 0.5)
    ############################################
    minimize_option = deepcopy(option)
    if n_chunk <= 1 :
        print "Input pdb < 150 residues, no slicing is required."
        print "Start minimizing the full pdb..."
        minimize_option.input_pdb = "before_min.pdb"
        minimize_option.out_pdb = "after_min.pdb"
        minimize_option.log_out = 'full_minimize_temp.out'
        erraser_minimize( minimize_option )
        print exists('after_min.pdb')
        move('after_min.pdb', option.out_pdb)
    else :
        print "Input pdb >= 150 residus, slice into %s chunks and minimize each one sequentially." % n_chunk
        res_slice_list = pdb_slice_into_chunks(option.input_pdb, n_chunk)
        current_chunk = 0
        for res_slice in res_slice_list :
            print "Start minimizing chunk %s..." % current_chunk
            print "Chunk %s residues: %s" % (current_chunk, res_slice)
            current_chunk += 1
            minimize_option.input_pdb = "before_min.pdb"
            minimize_option.out_pdb = "after_min.pdb"
            minimize_option.res_slice = res_slice
            minimize_option.log_out = "full_minimize_temp_%s.out" % current_chunk
            erraser_minimize( minimize_option )
            copy('after_min.pdb', 'before_min.pdb')
            print "Minimization for chunk %s ends sucessfully." % current_chunk
        move('after_min.pdb', option.out_pdb)

    os.chdir(base_dir)

    if not option.kept_temp_folder :
        remove(temp_dir)

    total_time=time.time()-start_time

    print '\n', "DONE!...Total time taken= %f seconds" % total_time
    print '###################################'
    if sys.stdout != sys.__stdout__:
        sys.stdout.close()
    if sys.stderr != sys.__stderr__:
        sys.stderr.close()
    sys.stdout = stdout
    sys.stderr = stderr
##### full_struct_slice_and_minimize end   ############################




##### seq_rebuild start ###############################################
def seq_rebuild( option ) :

    stdout = sys.stdout
    stderr = sys.stderr
    if option.log_out != "" :
        sys.stdout = open(option.log_out, 'w')
    if option.log_err != "" :
        sys.stderr = open(option.log_err, 'w')

    print '###################################'
    print 'Starting seq_rebuild...'
    start_time=time.time()
    option.finalize()

    #####Set temp folder#######################
    base_dir = os.getcwd()
    temp_dir = '%s/%s/' % (base_dir, basename(option.input_pdb).replace('.pdb', '_seq_rebuild_temp') )
    if exists(temp_dir) :
        print 'Temporary directory %s exists... Remove it and create a new folder.' % temp_dir
        remove(temp_dir)
        os.mkdir(temp_dir)
    else :
        print 'Create temporary directory %s...' % temp_dir
        os.mkdir(temp_dir)

    print '###################################'
    #####################################################
    os.chdir(temp_dir)

    copy(option.input_pdb, 'temp.pdb')

    sucessful_res = []
    failed_res = []
    SWA_option = deepcopy(option)
    for res in option.rebuild_res_list :
        print 'Starting to rebuild residue %s' % res
        SWA_option.input_pdb = 'temp.pdb'
        SWA_option.rebuild_res = res
        SWA_option.log_out = 'seq_rebuild_temp_%s.out' % res

        SWA_rebuild_erraser( SWA_option )

        print 'Job completed for residue %s' % res

        rebuilt_pdb_merge0 = './temp_pdb_res_%d/output_pdb/S_000000_merge.pdb' % res
        rebuilt_pdb_orig0 = './temp_pdb_res_%d/output_pdb/S_000000.pdb' % res
        rebuilt_pdb_merge1 = './temp_pdb_res_%d/output_pdb/S_000001_merge.pdb' % res
        rebuilt_pdb_orig1 = './temp_pdb_res_%d/output_pdb/S_000001.pdb' % res
        rebuilt_pdb_final = ''

        if exists(rebuilt_pdb_merge0) :
            rebuilt_pdb_final = rebuilt_pdb_merge0
        elif exists(rebuilt_pdb_orig0) :
            rebuilt_pdb_final = rebuilt_pdb_orig0
        elif exists(rebuilt_pdb_merge1) :
            rebuilt_pdb_final = rebuilt_pdb_merge1
        elif exists(rebuilt_pdb_orig1) :
            rebuilt_pdb_final = rebuilt_pdb_orig1

        if rebuilt_pdb_final != '' :
            print "Residue %d is sucessfully rebuilt!" % res
            sucessful_res.append(res)
            copy(rebuilt_pdb_final, 'temp.pdb')
        else :
            print "No suitable alternative structure can be sampled."
            print "Residue %d is not rebuilt!" % res
            failed_res.append(res)

        if not option.verbose :
            remove('temp_pdb_res_%d' % res)

        print '###################################'

    copy('temp.pdb', option.out_pdb)
    os.chdir(base_dir)

    if not option.kept_temp_folder :
        remove(temp_dir)

    print "All rebuilding moves completed sucessfully!!!!"
    print 'sucessful_res: %s' % sucessful_res
    print 'failed_res: %s' % failed_res

    total_time=time.time()-start_time
    print '\n', "DONE!...Total time taken= %f seconds" % total_time
    print '###################################'
    if sys.stdout != sys.__stdout__:
        sys.stdout.close()
    if sys.stderr != sys.__stderr__:
        sys.stderr.close()
    sys.stdout = stdout
    sys.stderr = stderr
##### seq_rebuild end   ###############################################




##### SWA_rebuild_erraser start #######################################
def SWA_rebuild_erraser( option ) :

    stdout = sys.stdout
    stderr = sys.stderr
    if option.log_out != "" :
        sys.stdout = open(option.log_out, 'w')
    if option.log_err != "" :
        sys.stderr = open(option.log_err, 'w')

    print '###################################'
    print 'Starting SWA_rebuild_erraser...'
    start_time=time.time()

    native_screen = True
    if option.native_screen_RMSD > 10.0 :
        native_screen = False

    ##############location of executable and database:###########
    database_folder = rosetta_database_path(option.rosetta_database)
    rna_swa_test_exe = rosetta_bin_path("swa_rna_main", option.rosetta_bin )
    #############folder_name###################################
    base_dir=os.getcwd()

    main_folder =          abspath("%s_res_%d" %(basename(option.input_pdb).replace('.','_'), option.rebuild_res) )
    temp_folder=           abspath("%s/temp_folder/" %(main_folder) )
    sampling_folder=       abspath("%s/sampling/" %(main_folder) )
    cluster_folder=        abspath("%s/cluster/" %(main_folder))
    output_pdb_folder=     abspath("%s/output_pdb/" %(main_folder))
    precluster_pdb_folder= abspath("%s/precluster_pdb/" %(main_folder))

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

    native_pdb = abspath("%s/native_struct.pdb" % output_pdb_folder)
    copy(option.input_pdb, native_pdb)
    #####################Slice out the rebuild region############
    cutpoint_final = []
    res_sliced_all = []
    native_pdb_final = native_pdb
    rebuild_res_final = option.rebuild_res
    if option.slice_nearby :
        native_pdb_final = native_pdb.replace('native_struct.pdb', 'native_struct_sliced.pdb')

        rebuilding_res = []
        rebuilding_res.append( option.rebuild_res )
        if option.rebuild_res != 1 and (not option.rebuild_res - 1 in option.cutpoint_rs) :
            rebuilding_res.append( option.rebuild_res - 1 )
        if option.rebuild_res != get_total_res(native_pdb) and (not option.rebuild_res in option.cutpoint_rs) :
            rebuilding_res.append( option.rebuild_res + 1 )

        res_sliced_all = pdb_slice_with_patching( native_pdb, native_pdb_final, rebuilding_res ) [0]

        rebuild_res_final = res_sliced_all.index(option.rebuild_res) + 1
        for cutpoint in option.cutpoint_rs :
            if cutpoint in res_sliced_all :
                cutpoint_final.append( res_sliced_all.index(cutpoint) + 1 )
    else :
        cutpoint_final = option.cutpoint_rs

    total_res = get_total_res(native_pdb_final)

    print "rebuilding_res = %s" % rebuild_res_final
    print "res_sliced = %s" % res_sliced_all
    print "cutpoint_final = %s" % cutpoint_final
    print "total_res = %d " % total_res

    #######Decide whether to sample syn pyrimidine################
    if option.allow_syn_pyrimidine and option.search_syn_pyrimidine_only_when_native_syn :
        chi_angle = find_chi_angle(native_pdb_final, rebuild_res_final)
        is_syn = ( chi_angle < 140 and chi_angle > -40 )
        if is_syn :
            allow_syn_pyrimidine = 'true'
        else :
            allow_syn_pyrimidine = 'false'
    else :
        allow_syn_pyrimidine =  str(option.allow_syn_pyrimidine).lower()
    ##########Check if the rebuilding Rsd is at chain break###########
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
        option.ideal_geometry = False

    start_pdb = ""
    if option.ideal_geometry :
        start_pdb = "%s/missing_rebuild_res_native.pdb" % temp_folder
        pdb_slice(native_pdb_final, start_pdb, "%d-%d %d-%d" % (1, rebuild_res_final-1, rebuild_res_final+1, total_res) )
    else :
        start_pdb = native_pdb_final

    #####################Create fasta file########################
    fasta_file=temp_folder + '/fasta'
    pdb2fasta(native_pdb_final, fasta_file, using_protein = option.rna_prot_erraser)
    #########################Common Options######################
    common_cmd = ""

    common_cmd += " -database %s " % database_folder

    common_cmd += " -VERBOSE %s" % str(option.verbose).lower()

    common_cmd += " -fasta %s " % fasta_file

    common_cmd += " -input_res "
    for n in range(1,total_res+1) :
        if not (option.ideal_geometry and n == rebuild_res_final) :
            common_cmd += "%d " %n

    common_cmd+= " -fixed_res "
    for n in range(1,total_res+1) :
        if n != rebuild_res_final :
            common_cmd += "%d " %n

    #PHENIX conference -- HACK -- try to specify exactly the jump points. Needed for RNA/protein poses.
    #protein case
    if option.rna_prot_erraser :
        common_cmd += " -jump_point_pairs %d-%d " % ( rebuild_res_final-1, rebuild_res_final+1 )
    else : #RNA only original case
        common_cmd += " -jump_point_pairs NOT_ASSERT_IN_FIXED_RES 1-%d " % total_res

    common_cmd += " -alignment_res 1-%d " % total_res

    common_cmd += " -rmsd_res %d " %(total_res)
    common_cmd += " -native " + native_pdb_final
    common_cmd += " -score:weights %s " % option.scoring_file
    
    #Rescue 2012 defaults 
    if option.enlarge_H_lj is True:
        common_cmd += " -enlarge_H_lj %s " % str(option.enlarge_H_lj).lower()
   if option.enlarge_lj_hbond_dis is True:
        common_cmd += " -lj_hbond_hdis 1.95 "
        common_cmd += " -lj_hbond_OH_donor_dis 3.0 "
    if option.restore_pre_talaris_2013_behavior is True:
        common_cmd += " -restore_pre_talaris_2013_behavior %s " % str(option.restore_pre_talaris_2013_behavior).lower()
        common_cmd += " -score::smooth_fa_elec %s " % str(option.smooth_fa_elec).lower()
    if option.use_2prime_OH_potential is False:
        common_cmd += " -use_2prime_OH_potential %s " % str(option.use_2prime_OH_potential).lower()

    if option.map_file != "" :
        common_cmd += " -edensity:mapfile %s " % option.map_file
        common_cmd += " -edensity:mapreso %s " % option.map_reso
        common_cmd += " -edensity:realign no "

    if len(cutpoint_final) != 0 :
        common_cmd += " -cutpoint_open "
        for cutpoint in cutpoint_final :
            common_cmd += '%d ' % cutpoint

    if option.fcc2012_new_torsional_potential :
        common_cmd += " -score:rna_torsion_potential FCC2012_RNA11_based_new "
    elif option.new_torsional_potential :
        common_cmd += " -score:rna_torsion_potential RNA11_based_new "

    common_cmd += " -rna::corrected_geo %s " % str(option.corrected_geo).lower()
    common_cmd += " -rna::rna_prot_erraser %s " % str(option.rna_prot_erraser).lower()
    ################Sampler Options##################################
    sampling_cmd = rna_swa_test_exe + ' -algorithm rna_sample '
    if not is_chain_break :
        sampling_cmd += '-erraser true '
        sampling_cmd += '-sampler_extra_epsilon_rotamer true '

    sampling_cmd += " -s %s " % start_pdb
    sampling_cmd += " -out:file:silent blah.out "
    sampling_cmd += " -output_virtual true "
    sampling_cmd += " -rm_virt_phosphate true "
    sampling_cmd += " -sampler_extra_chi_rotamer true "
    sampling_cmd += " -cluster::radius %s " % 0.3
    sampling_cmd += " -centroid_screen true "
    #sampling_cmd += " -VDW_atr_rep_screen false "
    sampling_cmd += " -sampler_allow_syn_pyrimidine %s " % allow_syn_pyrimidine
    sampling_cmd += " -minimize_and_score_native_pose %s " % str(option.include_native).lower()
    sampling_cmd += " -native_edensity_score_cutoff %s " % option.native_edensity_cutoff
    sampling_cmd += " -constraint_chi %s " % str(option.constrain_chi).lower()
    if native_screen:
        sampling_cmd += " -rmsd_screen %s " % option.native_screen_RMSD
    sampling_cmd += " -sampler_num_pose_kept %s " % option.num_pose_kept
    sampling_cmd += " -PBP_clustering_at_chain_closure true "
    sampling_cmd += " -allow_chain_boundary_jump_partner_right_at_fixed_BP true "
    sampling_cmd += " -allow_virtual_side_chains false"
    sampling_cmd += " -sampler_perform_phosphate_pack false"
    sampling_cmd += " -add_virt_root true "


    ##################################################################
    if not is_chain_break :
        #################Use Analytical Loop Closure##################
        if(exists(sampling_folder)):
            print "warning...sampling_folder:%s already exist...removing it...! " % sampling_folder
            remove(sampling_folder)
        os.mkdir(sampling_folder)

        if option.is_append :
            print  '\n', "Rebuilding res %s by attaching to res %s" % (rebuild_res_final, rebuild_res_final-1),  '\n'
        else :
            print  '\n', "Rebuilding res %s by attaching to res %s" % (rebuild_res_final, rebuild_res_final+1),  '\n'

        os.chdir( sampling_folder )

        specific_cmd =""
        specific_cmd += " -sample_res %d " % rebuild_res_final
        if option.is_append :
            specific_cmd += " -cutpoint_closed %d " % rebuild_res_final
        else :
            specific_cmd += " -cutpoint_closed %d " % (rebuild_res_final - 1)

        command = sampling_cmd + ' ' + specific_cmd + ' ' + common_cmd
        if (option.verbose): print  '\n', command, '\n'

        subprocess_call( command, 'sampling_1.out', 'sampling_1.err' )

        os.chdir( base_dir )

    #########Sample chainbreak residues with original SWA rebuild#####
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

        if(option.verbose): print  '\n', command, '\n'

        subprocess_call( command, 'sampling_1.out', 'sampling_1.err' )


        os.chdir( base_dir)

    ###################Clustering############
    #Just output the lowest energy decoy instead of clustering if num_pose_kept_cluster = 1
    is_clustering = (option.num_pose_kept_cluster != 1)
    if is_clustering :
        print '\n',"ALMOST DONE...sorting/clustering/extracting output_pdb....", '\n'
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
        cluster_args += " -in:file:silent_struct_type rna"
        cluster_args += " -in:file:silent %s/blah.out " % sampling_folder
        cluster_args += " -PBP_clustering_at_chain_closure true "
        cluster_args += " -allow_chain_boundary_jump_partner_right_at_fixed_BP true "

        no_clustering  = " -suite_cluster_radius 0.0 "
        no_clustering += " -loop_cluster_radius 0.0 "

        if(option.verbose):  ##This is just for control purposes...
            command = (cluster_args + ' ' + common_cmd + no_clustering +
                " -recreate_silent_struct false  -out:file:silent %s" % CONTROL_filename)

            if exists(sampling_folder + '/blah.out'):
                print '\n', command ,'\n'
                subprocess_call( command, 'CONTROL.out', 'CONTROL.err' )

            #This one is with alignment to native_pose...however wary that SCORE output is not correct...
            command = (cluster_args + ' ' + common_cmd +  no_clustering +
              " -recreate_silent_struct true -out:file:silent %s" % precluster_filename)

            if exists(sampling_folder + '/blah.out'):
                if(option.verbose): print '\n', command ,'\n'
                subprocess_call( command, 'precluster.out', 'precluster.err' )

        with_clustering  = ""
        with_clustering += " -suite_cluster_radius %s " % option.cluster_RMSD
        with_clustering += " -loop_cluster_radius 999.99 "
        with_clustering += " -clusterer_num_pose_kept %d " % option.num_pose_kept_cluster

        command = (cluster_args + ' ' + common_cmd +  with_clustering +
          " -recreate_silent_struct true -out:file:silent %s" % cluster_filename )

        if exists(sampling_folder + '/blah.out'):
            if(option.verbose): print '\n', command ,'\n'
            subprocess_call( command, 'cluster.out', 'cluster.err' )

        os.chdir( base_dir)

        if option.verbose :
            if exists(CONTROL_filename) :
                score_line = subprocess_out('head -n 2 %s ' % CONTROL_filename ) [1]
                subprocess_call('echo "%s"' % score_line, '%s/output_pdb_CONTROL.txt' % main_folder)
                subprocess_call("grep SCORE %s | sort -nk2" % CONTROL_filename, "%s/output_pdb_CONTROL.txt" % main_folder)

            if exists(precluster_filename) :
                score_line = subprocess_out('head -n 2 %s ' % precluster_filename ) [1]
                subprocess_call('echo "%s"' % score_line, '%s/output_pdb_precluster.txt' % main_folder)
                subprocess_call("grep SCORE %s | sort -nk2" % precluster_filename, "%s/output_pdb_precluster.txt" % main_folder)

        if exists(cluster_filename) :
            score_line = subprocess_out('head -n 2 %s ' % cluster_filename ) [1]
            subprocess_call('echo "%s"' % score_line, '%s/output_pdb.txt' % main_folder)
            subprocess_call("grep SCORE %s | sort -nk2" % cluster_filename, "%s/output_pdb.txt" % main_folder)
        else :
            output = open("%s/output_pdb.txt" % main_folder, 'w')
            output.write("No silent file is being output during rebuilding")

        if option.verbose and exists(precluster_filename) :
            extract_pdb(precluster_filename, precluster_pdb_folder, option.rosetta_bin, option.rosetta_database, rna_prot_erraser=option.rna_prot_erraser)
        if exists(cluster_filename) :
            extract_pdb(cluster_filename, output_pdb_folder, option.rosetta_bin, option.rosetta_database, rna_prot_erraser=option.rna_prot_erraser)

    else : #Fast outputing
        out_slient_file = "%s/blah.out"  % sampling_folder
        if exists(out_slient_file) :
            extract_pdb(out_slient_file, output_pdb_folder, option.rosetta_bin, option.rosetta_database, extract_first_only = True, rna_prot_erraser=option.rna_prot_erraser)

    ##############Merge the sliced region back to starting pdb######
    if option.slice_nearby :
        os.chdir( output_pdb_folder )
        pdb_file_list = glob("S*.pdb")
        for pdb_file in pdb_file_list :
            sliced2orig_merge_back( native_pdb, pdb_file, pdb_file.replace('.pdb', '_merge.pdb'), res_sliced_all )
    #############################################################
    os.chdir(base_dir)

    if not option.verbose :
        remove(temp_folder)
        remove(sampling_folder)
        if is_clustering :
            remove(cluster_folder)
            remove(precluster_pdb_folder)
            if exists(cluster_filename):
                remove(cluster_filename)

    total_time=time.time()-start_time

    print '\n', "DONE!...Total time taken= %f seconds" % total_time , '\n'
    print '###################################'
    if sys.stdout != sys.__stdout__:
        sys.stdout.close()
    if sys.stderr != sys.__stderr__:
        sys.stderr.close()
    sys.stdout = stdout
    sys.stderr = stderr
##### SWA_rebuild_erraser end   #######################################
