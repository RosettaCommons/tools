from os.path import exists, expandvars, abspath
from math import sqrt
import subprocess
import sys
import os
import string
import os.path

#####################################################
def rosetta_bin(exe_file = "") :
    """
    Return the absolute path to the given exexutable file in Rosetta.
    If no input is given, return the Rosetta bin folder path.
    """
    rosetta_folder = expandvars("$ROSETTA")
    if rosetta_folder == "$ROSETTA" :
        error_exit_with_message("USER need to set environmental variable $ROSETTA and pointed it to the Rosetta folder!")

    exe_folder = rosetta_folder + "/rosetta_source/bin/"
    check_path_exist(exe_folder)
    exe_path = exe_folder + exe_file
    check_path_exist(exe_path)
    return exe_path
#####################################################
def rosetta_database() :
    """
    Return the absolute path to the Rosetta database.
    If no input is given, return the Rosetta bin folder path.
    """
    rosetta_folder = expandvars("$ROSETTA")
    if rosetta_folder == "$ROSETTA" :
        error_exit_with_message("USER need to set environmental variable $ROSETTA and pointed it to the Rosetta folder!")

    database_folder = rosetta_folder + "/rosetta_database/"
    check_path_exist(database_folder)
    return database_folder    
#####################################################
def error_exit_with_message(message):
    print >> sys.stdout, "A error had occur, see errfile.txt"
    print >> sys.stderr, "#####################################################################"
    print >> sys.stderr, "Error!: "  + message
    print >> sys.stderr, "#####################################################################"

    sys.stdout.flush()
    sys.stderr.flush()
    assert(False)
###################################################
def subprocess_call(command):
    """
    Execute the given commands in /usr/bin/sh.
    Error exit if failed.
    """
    sys.stdout.flush()
    sys.stderr.flush()

    subprocess.check_call(command, shell=True)

    sys.stdout.flush()
    sys.stderr.flush()
###################################################
def subprocess_out(command):
    """
    Execute the given commands in /usr/bin/sh.
    Error exit if failed.
    Return the output as a list of each lines.
    """
    sys.stdout.flush()
    sys.stderr.flush()
    
    out_list = subprocess.check_output(command, shell=True).split('\n')

    sys.stdout.flush()
    sys.stderr.flush()
    return out_list
###################################################
def get_total_res(pdbname):
    """
    Return the total number of residues ina pdb file.
    """
    netpdbname = pdbname
    check_path_exist(netpdbname)

    oldresnum = '   '
    count = 0;
    for line in open(netpdbname):
        if len(line) > 4 and line[0:4] == 'ATOM':
            resnum = line[22:26]
            if resnum != oldresnum :
                count = count + 1
                oldresnum = resnum

    return count
################################################################
def parse_options( argv, tag, default):
    """
    Parse various input options.
    """
    tag_argv = '-%s' % tag
    if tag_argv in argv :
        pos = argv.index( tag_argv )   ###Position of the option name 

        if ( ( pos == (len( argv )-1) or argv[ pos+1 ][0] == '-' ) ):
            error_exit_with_message("Invalid parse_option input")

        if (default == "False" or default == "True"): # Python boolean
            if ( argv[ pos + 1 ] == "True"):
                value=True
            elif ( argv[ pos + 1 ] == "False"):
                value=False
            else: 
                error_exit_with_message('(%s != "True") and  (%s != "False")' % (tag, tag)) 
        elif( isinstance( default, str ) ): 
            value = argv[ pos + 1 ]
        elif isinstance( default, int ):
            value = int( argv[ pos + 1 ] )
        elif isinstance( default, float ):
            value = float( argv[ pos + 1 ] )
        else:
            error_exit_with_message("Invalid parse_option default")
            
        print "%s = %s" %(tag, argv[ pos + 1 ])            
        return value    
            
    else: #Return the default value
        print "%s = %s" %(tag, default)
        if(default=="False" or default=="True"): 
            if( default == "True"):
                return True
            elif( default == "False"):
                return False
            else: 
                error_exit_with_message('(%s != "True") and  (%s != "False")' % (tag, tag)) 
        else:
            return default
###############################
def parse_option_int_list ( argv, tag ) :
    """
    Parse input cmdline option with the following format:
    ... -tag 1 5 7-20 40 45-60 ...
    Return a list of input numbers
    """
    list_load = []
    tag_argv = '-%s' % tag

    if tag_argv in argv :
        pos = argv.index(tag_argv)

        for i in range( pos + 1, len(argv) ) :
            input_string = argv[i]
            if input_string[0] == '-' :
                break
            split_string = input_string.split('-')

            for num in split_string :
                try :
                    int(num)
                except :
                    error_exit_with_message("Incorrect input for -%s: instance %s" % (tag, input_string) )

            if len(split_string) == 1 :
                list_load.append(int(split_string[0]))
            elif len(split_string) == 2 :
                start_num = int(split_string[0])
                end_num = int(split_string[1])
                for j in range(start_num, end_num + 1) :
                    list_load.append(j)
            else :
                error_exit_with_message("Incorrect input for -%s: instance %s" % (tag, input_string) )
    print "%s = %s" % (tag, list_load)        
    return list_load
###############################
def parse_option_string_list ( argv, tag ) :
    """
    Parse input cmdline option with the following format:
    ... -tag ss bbb fef ...
    Return a list of input strings
    """
    list_load = []
    tag_argv = '-%s' % tag

    if tag_argv in argv :
        pos = argv.index(tag_argv)

        for i in range( pos + 1, len(argv) ) :
            input_string = argv[i]
            if input_string[0] == '-' :
                break
            list_load.append(input_string)
    print "%s = %s" % (tag, list_load)        
    return list_load
##########################################
def parse_option_chain_res_list ( argv, tag ) :
    """
    Parse input cmdline option with the following format:
    ... -tag A10 B20-22...
    Return a list of input strings with number expanded:
    ['A10', 'B20', 'B21', 'B22', ... ]
    """
    list_load = []
    tag_argv = '-%s' % tag
    
    if tag_argv in argv :
        pos = argv.index(tag_argv)

        for i in range( pos + 1, len(argv) ) :
            input_string = argv[i]
            if input_string[0] == '-' :
                break
            input_string = input_string.upper()
            if not input_string[0] in string.uppercase :
                error_exit_with_message( "Incorrect input for -%s: instance %s" % (tag, input_string) )
            
            chain_id = input_string[0]
            input_string = input_string[1:]
            split_string = input_string.split('-')

            for num in split_string :
                try :
                    int(num)
                except :
                    error_exit_with_message("Incorrect input for -%s: instance %s" % (tag, argv[i]) )

            if len(split_string) == 1 :
                list_load.append('%s%s' % (chain_id, split_string[0]) )
            elif len(split_string) == 2 :
                start_num = int(split_string[0])
                end_num = int(split_string[1])
                for j in range(start_num, end_num + 1) :
                    list_load.append('%s%s' % (chain_id, j))
            else :
                error_exit_with_message("Incorrect input for -%s: instance %s" % (tag, argv[i]) )
    print "%s = %s" % (tag, list_load)        
    return list_load
#####################################################
def rosetta_ready_set( input_pdb, out_name ) : 
    """
    Call Rosetta to read in a pdb and output the model right away.
    Can be used to ensure the file have the rosetta format for atom name, ordering and phosphate OP1/OP2.
    """
    check_path_exist(input_pdb)

    command = rosetta_bin("erraser_minimizer.linuxgccrelease")
    command += " -database %s" % rosetta_database()
    command += " -native %s" % input_pdb
    command += " -out_pdb %s" % out_name
    command += " -ready_set_only True"  
    command += " > ready_set.out"  
    subprocess_call( command )
    return True
#####################################################
def extract_pdb(silent_file, output_folder_name):
    """
    Extract pdb's from Rosetta silent files.
    """

    check_path_exist( silent_file )

    remove_variant_types = "false"

    #find the first discription line 

    ##########################
    line=""

    infile=open(silent_file,"r")

    found_description_line=False

    while (found_description_line == False):


        line=infile.readline()
        if line.find('desc') != -1: #found the description line
            found_description_line=True  

    line_list=line.split()

    infile.close();

    ########################

    column_num=0

    while column_num < len(line_list):
        column_num += 1
        if  line_list[column_num].find('desc') != -1:
            break
    #########################################
    tags_string = ""

    for line in open(silent_file):

        if line.find('SCORE') == -1: #Line contain the word SCORE
            continue

        if line.find('desc') != -1: #Line is not a column_name line
            continue

        tag=line.split()[column_num]

        tags_string += ' ' + tag

    #########################################################

    base_dir = os.getcwd()

    if not exists( output_folder_name ) :
        os.mkdir(output_folder_name)
    
    os.chdir( output_folder_name )

    command = rosetta_bin("rna_extract.linuxgccrelease")

    command += " -tags " + tags_string
    command += " -in:file:silent " + silent_file 
    command += " -in:file:silent_struct_type  binary_rna"
    command += " -database %s" % rosetta_database()
    command += " -remove_variant_cutpoint_atoms " + remove_variant_types 


    command += ' > extraction_output.txt'
    subprocess_call( command )

    os.chdir( base_dir )
    return True
#####################################################
def find_error_res(input_pdb) :
    """
    Return a list of error resdiue in the given pdb file (RNA only).
    Use phenix.rna_validate.
    """
    check_path_exist(input_pdb)
    output = ""
    try :
        output = subprocess_out("phenix.rna_validate %s" % input_pdb)
    except :
        error_exit_with_message("Error in find_error_res!!!!")

    error_types = ["Pucker", "Bond", "Angle", "Suite"]
    current_error = 0
    line = -1
    error_res = [ [], [], [], [] ] #Pucker, Bond, Angle, Suite error_res

    while current_error != 4 or line < len(output) - 1 :
        line += 1
        if len( output[line] ) < 7 :
            continue
        if error_types[current_error] in output[line] :
            line += 2
            while len( output[line] ) > 7 :
                res = int( output[line].split() [2]) 
                error_res[current_error].append( res )
                line += 1
            current_error += 1

    suite_res = []
    for res in error_res[3] :
        if res - 1 > 0 :
           suite_res.append(res - 1)
    error_res.append(suite_res)

    error_res_final = []
    for res_list in error_res :
        for res in res_list :
            if not res in error_res_final :
                error_res_final.append(res)
    error_res_final.sort()
    return error_res_final
#############################################
def pdb2fasta(input_pdb, fasta_out) :
    """
    Extract fasta info from pdb.
    """
    check_path_exist(input_pdb)

    longer_names={' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u',
                  'ADE': 'a', 'CYT': 'c', 'GUA': 'g', 'URA': 'u'}
    output = open(fasta_out, 'w')
    output.write( ">%s\n" % os.path.basename(input_pdb) )
    oldresnum = ''
    for line in open(input_pdb) :
        if len(line) > 20 : 
            if line[0:3] == 'TER':
                output.write('\n')

            if line[0:4] == 'ATOM':
                resnum = line[23:26]
                if not resnum == oldresnum:
                    longname = line[17:20]
                    if longer_names.has_key(longname):
                        output.write( longer_names[longname] );
                    else:
                        output.write( 'X')
                oldresnum = resnum
    output.write('\n')
    output.close()
    return True
##################################################
def pdb_slice(input_pdb, out_name, segment) :
    """
    Slice the input pdb file with given segment.
    Segment format example: '1-3 8 10-30 45 46 40'
    Or input a python int list as the segment.
    Return list of cut residue
    """
    check_path_exist(input_pdb)
    
    kept_res = []
    if type(segment) is list :
        kept_res = segment
    elif type(segment) is str :
        for elem in segment.split() :
            if '-' in elem :
                res_range = elem.split('-')
                for i in range(int(res_range[0]), int(res_range[1]) + 1) :
                    kept_res.append(i)
            else :
                kept_res.append( int(elem) )
    
    if len(kept_res) == 0 :
        error_exit_with_message("Invalid input of 'segment' option!")


    output = open(out_name, 'w')
    previous_res = -1
    old_res = 0
    new_res = 0
    new_atom = 0
    for line in open(input_pdb):
        if len(line) > 20 and line[0:4] == 'ATOM' :
            current_res = int(line[22:26])
            if current_res != previous_res :
                old_res += 1
                previous_res = current_res
                if old_res in kept_res :
                    new_res += 1

            if old_res in kept_res :
                new_atom += 1
                output.write('%s%7d%s%4d%s' % (line[0:4], new_atom, line[11:22], new_res, line[26:]) )

    output.close()
    return kept_res
#####################################
def check_path_exist(path_name) :
    if exists(path_name) :
        return True
    else :
        error_exit_with_message("Path %s does not exist!" % path_name)
        return False
#####################################
def compute_squared_dist(coord1, coord2) :
    """
    compute the squared distance between two xyz vector (3D).
    """
    sq_dist = 0
    sq_dist += (coord1[0] - coord2[0]) * (coord1[0] - coord2[0])
    sq_dist += (coord1[1] - coord2[1]) * (coord1[1] - coord2[1])
    sq_dist += (coord1[2] - coord2[2]) * (coord1[2] - coord2[2])
    return sq_dist
####################################
def compute_dist(coord1, coord2) :
    """
    compute the distance between two xyz vector (3D).
    """
    sq_dist = compute_squared_dist(coord1, coord2)
    return sqrt(sq_dist)
####################################
def load_pdb_coord(input_pdb) :
    """
    Load in the pdb and return the coordinates for each heavy atom in each residue.
    Also return a list of C1* coordinates
    """
    check_path_exist(input_pdb)

    current_res = 1
    coord_all = []
    coord_res = []
    coord_C1 = []
    for line in open(input_pdb) :
        if len(line) < 22 :
            continue
        if line[0:4] != 'ATOM' :
            continue
        if line[13] == 'H' or line[12] == 'H' or line[77] == 'H':
            continue
        res = int(line[22:26])
        if res != current_res :
            coord_all.append(coord_res)
            current_res = res
            coord_res = []

        coord_cur = [float(line[31:38]), float(line[39:46]), float(line[47:54])]
        coord_res.append(coord_cur)
        if line[13:16] == 'C1*' :
            coord_C1.append( coord_cur )
    coord_all.append(coord_res)
    if len(coord_C1) != len(coord_all) : 
        error_exit_with_message("Number of residues != number of C1*!!!")

    return [coord_all, coord_C1]
####################################
def find_nearby_res(input_pdb, input_res, dist_cutoff, is_reload = True):
    """
    Find nearby residues to the input residues by a distance_cutoff.
    All the residues should be in the same chain with continous numbering starting from 1.
    """
    check_path_exist(input_pdb)

    try :
        coord_all
        coord_C1
        if is_reload :
            [coord_all, coord_C1] = load_pdb_coord(input_pdb)
    except :
        [coord_all, coord_C1] = load_pdb_coord(input_pdb)
        
    res_list = []
    for i in input_res :
        if not i in range(1, len(coord_all) + 1) :
            error_exit_with_message("Input redidues outside the range of pdb residues!")
        for j in range(1, len(coord_all) + 1) :
            if (j in input_res or j in res_list) : continue
            dist_C1 = compute_dist( coord_C1[i-1], coord_C1[j-1] )
            if dist_C1 > dist_cutoff + 8 :
                continue
            for coord_target_atom in coord_all[i-1] :
                is_break = False
                for coord_all_atom in coord_all[j-1] :
                    dist = compute_dist( coord_target_atom, coord_all_atom)
                    if dist < dist_cutoff :
                        res_list.append(j)
                        is_break = True
                        break
                if is_break :
                    break

    res_list.sort()
    return res_list
#############################################################
def rosetta2phenix_merge_back(orig_pdb, rs_pdb, out_name) :
    """
    Merge the rosetta pdb back to the phenix refined pdb.
    """
    if orig_pdb == "" or rs_pdb == "" or out_name == "" :
        error_exit_with_message("Input argument imcomplete!")
    check_path_exist(orig_pdb)
    check_path_exist(rs_pdb)

    def is_line_rna (line) :
        res_name = line[17:20];
        if (res_name == "  G") : return True;
        if (res_name == "  C") : return True;
        if (res_name == "  A") : return True;
        if (res_name == "  U") : return True;

        if (res_name == " rG") : return True;
        if (res_name == " rC") : return True;
        if (res_name == " rA") : return True;
        if (res_name == " rU") : return True;

        return False;

    ###################################
    output_pdb = open(out_name, 'w')
    lines_rs = open(rs_pdb).readlines()
    lines_orig = open(orig_pdb).readlines()
    atom_name_convert = { " O1P":" OP1", " O2P":" OP2", " O3P":" OP3", 
                          "1H5'":" H5'", "2H5'":"H5''", "1H2'":" H2'", 
                          "1H2 ":" H21", "2H2 ":" H22", "2HO'":"HO2'",
                          "3HO'":"HO3'", "5HO'":"HO5'", "1H4 ":" H41", 
                          "2H4 ":" H42", "1H6 ":" H61", "2H6 ":" H62" }

    current_res = ''
    current_chain = ''
    res_rs = ''
    atom_index = 0
    for line_orig in lines_orig :
        header = line_orig[0:6] 
        if header == 'ATOM  ' and is_line_rna(line_orig) :
            chain_orig = line_orig[21]
            res_name_orig = line_orig[17:20]
            res_orig = line_orig[22:27]
            atom_name_orig = line_orig[12:16]
            if (res_orig != current_res or chain_orig != current_chain) :
                lines_rs_temp = lines_rs[:]
                for line_rs in lines_rs_temp :
                    if (line_rs[0:4] != 'ATOM') : 
                        lines_rs.remove(line_rs)
                        continue
                    res = line_rs[22:27]
                    if res_rs != res :
                        res_rs = res
                        break
                    else :
                        atom_index += 1
                        atom_name = line_rs[12:16].replace('*', "'")
                        if atom_name_convert.has_key(atom_name) :
                            atom_name = atom_name_convert[atom_name]   
                        output_line = (line_rs[0:6] + str(atom_index).rjust(5) + ' ' + atom_name +
                        line_rs[16:21] + current_chain + current_res + line_rs[27:])
                        output_line = output_line.replace('rG',' G')
                        output_line = output_line.replace('rA',' A')
                        output_line = output_line.replace('rU',' U')
                        output_line = output_line.replace('rC',' C')
                        output_pdb.write(output_line)
                        lines_rs.remove(line_rs)
                current_res = res_orig
                current_chain = chain_orig

            atom_index += 1
            output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:])

            lines_rs_temp = lines_rs[:]
            for line_rs in lines_rs_temp :
                if (line_rs[0:4] != 'ATOM') : continue
                res = line_rs[22:27]
                res_name = line_rs[19].rjust(3)

                atom_name = line_rs[12:16].replace('*', "'")
                if atom_name_convert.has_key(atom_name) :
                    atom_name = atom_name_convert[atom_name]

                if (res_name != res_name_orig) : continue
                if (atom_name != atom_name_orig) : continue
                if (res != res_rs) : break
                output_line = (output_line[0:27] + line_rs[27:55] + output_line[55:])
                lines_rs.remove(line_rs)
                break

            output_pdb.write(output_line)

        elif header == 'ANISOU' : 
            output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:])
            output_pdb.write(output_line)
        elif (header == 'HETATM' or 
            (header == 'ATOM  ' and (not is_line_rna(line_orig) ) ) ):
            res_name_orig = line_orig[17:20]
            atom_index += 1
            output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:])
            output_pdb.write(output_line)
        elif (header != 'MASTER' and header[0:3] != 'END' and header != 'ANISOU') :
            output_pdb.write(line_orig)

    output_pdb.write('END                                ')
    return True

###############################################
def rosetta2std_pdb (input_pdb, out_name, cryst_line = "") :
    """
    Convert the rosetta pdb file to standard pdb format. (RNA only)
    Prepend the CRYST1 line if given in input.
    """
    check_path_exist(input_pdb)

    atom_name_convert = { " O1P":" OP1", " O2P":" OP2", " O3P":" OP3",
                          "1H5'":" H5'", "2H5'":"H5''", "1H2'":" H2'", 
                          "1H2 ":" H21", "2H2 ":" H22", "2HO'":"HO2'",
                          "3HO'":"HO3'", "5HO'":"HO5'", "1H4 ":" H41", 
                          "2H4 ":" H42", "1H6 ":" H61", "2H6 ":" H62" }

    res_name_convert = {" rG":"  G", " rA":"  A", " rU":"  U", " rC":"  C"}
    output = open(out_name, 'w')
    if cryst_line != "" :
        output.write("%s\n" % cryst_line)
    for line in open(input_pdb) :
        if len(line) <= 4 :
            continue
        elif line[0:4] != 'ATOM' :
            output.write(line)
            continue

        new_line = line.replace("*", "'")
             
        atom_name = new_line[12:16]
        if atom_name_convert.has_key(atom_name) :
            atom_name = atom_name_convert[atom_name]

        res_name = line[17:20]
        if res_name_convert.has_key(res_name) :
            res_name = res_name_convert[res_name]

        atom_type = line[13]
        if line[12] == 'H' :
            atom_type = 'H'

        new_line_list = list(new_line)
        new_line_list[12:16] = atom_name
        new_line_list[17:20] = res_name
        new_line_list[77] = atom_type

        output_line = string.join(new_line_list, '')
        output.write(output_line)

    output.close()
    return True
####################################
def pdb2rosetta (input_pdb, out_name, alter_conform = 'A', PO_dist_cutoff = 2.0) :
    """
    Convert regular pdb file to rosetta format.
    Return a list for converting original residues # into new ones,
    a list of residue convalently bonded to removed HETATM do they can be fixed in further refinement,
    a list of cutpoints (terminal residues) in the structure,
    and the CRYST1 line in the model.
    """

    check_path_exist(input_pdb)

    atom_name_convert = { " OP1":" O1P", " OP2":" O2P", " OP3":" O3P", 
                          " H5'":"1H5'", "H5''":"2H5'", " H2'":"1H2'", 
                          " H21":"1H2 ", " H22":"2H2 ", "HO2'":"2HO'",
                          "HO3'":"3HO'", "HO5'":"5HO'", " H41":"1H4 ", 
                          " H42":"2H4 ", " H61":"1H6 ", " H62":"2H6 " }

    res_name_convert = { "  G":" rG", "G  ":" rG", "GUA":" rG",
                         "  A":" rA", "A  ":" rA", "ADE":" rA",
                         "  U":" rU", "U  ":" rU", "URA":" rU",
                         "  C":" rC", "C  ":" rC", "CYT":" rC" }

    res_conversion_list = []
    fixed_res_list = []
    cutpoint_list = []
    CRYST1_line = ''
    O3prime_coord_pre = []
    O3prime_coord_cur = []
    P_coord_cur = []
    is_previous_het = False
    is_current_het = False
    previous_res = ""
    res_no = 0
    atm_no = 0

    output = open(out_name, 'w')
    for line in open(input_pdb) :
        if len(line) <= 21 :  
            continue

        if line[0:6] == 'CRYST1' :
            CRYST1_line = line[:-1]
        elif line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM' :
            res_name = line[17:20]
            if line[0:6] == 'ATOM  ' :
                if res_name_convert.has_key(res_name) :
                    res_name = res_name_convert[res_name]
                else :
                    continue

            current_res = line[21:26]

            if current_res != previous_res :
                previous_res = current_res
                is_previous_het = is_current_het
                if (not is_previous_het) and (not is_current_het) :
                    if len(O3prime_coord_cur) == 0 or len(P_coord_cur) == 0 :
                        if res_no - 1 > 0 :
                            cutpoint_list.append(res_no - 1)

                is_current_het = ( line[0:6] == 'HETATM' )
                if len(O3prime_coord_cur) != 0 and len(O3prime_coord_cur) != 3 :
                    error_exit_with_message("More than one O3' in a residue!!")
                else :
                    O3prime_coord_pre = O3prime_coord_cur
                O3prime_coord_cur = []
                P_coord_cur = []
                if not is_current_het :
                    res_no += 1
                    orig_res = '%s%s' % (line[21], line[22:27])
                    orig_res = orig_res.replace(' ', '')
                    res_conversion_list.append(orig_res)

            if line[16] != ' ' and line[16] != alter_conform :
                continue

            atom_name = line[12:16]

            if atom_name_convert.has_key(atom_name) :
                atom_name = atom_name_convert[atom_name]
            atom_name = atom_name.replace("'", "*")

            if atom_name == " O3*" :
                coord = line[30:54].split()
                for elem in coord :
                    O3prime_coord_cur.append(float(elem))

            if atom_name == " P  " :
                coord = line[30:54].split()
                for elem in coord :
                        P_coord_cur.append(float(elem))

                if (not is_current_het) and is_previous_het :
                    if res_no - 1 > 0 :
                        cutpoint_list.append( res_no - 1 )
                
                if len(O3prime_coord_pre) == 3 :
                    PO_dist = compute_dist( P_coord_cur, O3prime_coord_pre )
                    if PO_dist > PO_dist_cutoff :
                        if (not is_previous_het) and (not is_current_het) :
                            if res_no - 1 > 0 :
                                cutpoint_list.append( res_no - 1 )
                    else :
                        if is_current_het and ( not is_previous_het ) :
                            if res_no > 0 :
                                fixed_res_list.append( res_no )

                        if ( not is_current_het ) and is_previous_het :
                            if res_no > 0 :
                                fixed_res_list.append( res_no )

            if not is_current_het :
                atm_no += 1
                line_list = list(line)
                line_list[4:11] = str(atm_no).rjust(7)
                line_list[12:16] = atom_name
                line_list[16] = ' '
                line_list[17:20] = res_name
                line_list[21] = 'A'
                line_list[22:26] = str(res_no).rjust(4)
                line_list[55:60] = " 1.00"
                output.write( string.join(line_list, '') )

    output.close()
    return [res_conversion_list, fixed_res_list, cutpoint_list, CRYST1_line]
##################################################
def res_wise_rmsd(pdb1, pdb2) :
    """
    Compute the rmsd of each residue from two pdb file. No extra alignment is performed.
    Return a list of rmsd of each residue.
    Assuming the pdb files have the same ordering of res name / atom name.
    """
    check_path_exist(pdb1)
    check_path_exist(pdb2)

    coord_pdb1 = load_pdb_coord(pdb1) [0]
    coord_pdb2 = load_pdb_coord(pdb2) [0]

    if len(coord_pdb1) != len(coord_pdb2) :
        error_exit_with_message("Two pdbs have different # of residues!!!")
        
    res_rmsd_list = []
    for i in range(0, len(coord_pdb1)) :
        if len(coord_pdb1[i]) != len(coord_pdb2[i]) :
            error_exit_with_message("Residue %s have different # of atoms in the two pdbs!!!" % (i+1) )
        res_rmsd = 0.0
        for j in range(0, len(coord_pdb1[i])) :
            res_rmsd += compute_squared_dist(coord_pdb1[i][j], coord_pdb2[i][j])
        res_rmsd = sqrt(res_rmsd)
        res_rmsd_list.append([i+1, res_rmsd])

    return res_rmsd_list

##################################################
def pdb_slice_into_chunks(input_pdb, n_chunk = -1) :
    """
    Slice the pdb into smaller chunks.
    Return a list of res_list for each chunk.
    """

    if n_chunk == -1 :
        error_exit_with_message("Input argument imcomplete!")
    elif n_chunk == 1 :
        return []
    check_path_exist(input_pdb)

    def fill_gaps_and_remove_isolated_res(res_list, total_res) :
        res_list.sort()

        res_list.sort()
        res_remove = []
        for i in range(0, len(res_list)) :
            if i == 0 and res_list[i+1] - res_list[i] != 1 :
                res_remove.append(res_list[i])
            elif i == len(res_list) - 1 and res_list[i] - res_list[i-1] != 1 :
                res_remove.append(res_list[i])
            elif res_list[i] - res_list[i-1] != 1 and res_list[i+1] - res_list[i] != 1 :
                res_remove.append(res_list[i])

        for res in res_remove :
            res_list.remove(res)

        new_res = []
        for i in range(0, len(res_list)) :
            if i == len(res_list) - 1 :
                if total_res - res_list[i] > 0 and  total_res - res_list[i] <= 4 :
                    new_res += range(res_list[i] + 1, total_res + 1)
            elif res_list[i+1] - res_list[i] > 1 and  res_list[i+1] - res_list[i] <= 4 :
                new_res += range(res_list[i] + 1, res_list[i+1])
        res_list += new_res

        return res_remove
     

    total_res = get_total_res(input_pdb)
    res_list_sliced = []
    sliced_list_final = []
    res_list_unsliced = range(1, total_res + 1)
    res_list_current = []
    chunk_size = 0
    n_chunk_left = n_chunk
    while len(res_list_unsliced) != 0 :
        if len(res_list_current) == 0 :
            res_list_current.append( res_list_unsliced.pop(0) )
            chunk_size = int( len(res_list_unsliced) / n_chunk_left * 1.1)
            n_chunk_left -= 1
        
        res_list_new = find_nearby_res(input_pdb, res_list_current, 4.0, False)
        res_remove = []
        for res in res_list_new :
            if res in res_list_sliced :
                res_remove.append(res)
        for res in res_remove :
            res_list_new.remove(res)

        if len(res_list_new) == 0 and (not len(res_list_current) >= chunk_size * 0.7) :
            while True :
                res = res_list_unsliced.pop(0)
                if not res in res_list_current :
                    res_list_new.append(res)
                    break

        res_list_current += res_list_new

        if len(res_list_current) >= chunk_size or len(res_list_new) == 0 :
            res_list_current.sort()
            fill_gaps_and_remove_isolated_res(res_list_current, total_res)

            res_remove = []
            for res in res_list_current :
                if res in res_list_unsliced :
                    res_remove.append(res)
            for res in res_remove :
                res_list_unsliced.remove(res)
                    
            sliced_list_final.append( res_list_current )
            res_list_sliced += res_list_current
            if n_chunk_left == 1 :
                res_list_current = res_list_unsliced
                removed_res = fill_gaps_and_remove_isolated_res(res_list_current, total_res) 
                sliced_list_final.append( res_list_current )
                while len(removed_res) != 0 :
                    res_remove = []
                    for res in removed_res :
                        res_list_near = find_nearby_res(input_pdb, [res], 2.0, False)
                        for res_list in sliced_list_final :
                            is_break = False
                            for res_near in res_list_near :
                                if res_near in res_list:
                                    res_list.append(res)
                                    res_list.sort()
                                    res_remove.append(res)
                                    is_break = True
                                    break
                            if is_break :
                                break

                    for res in res_remove :
                        removed_res.remove(res)
                break
            else :
                res_list_current = []
                res_list_current.append( res_list_unsliced.pop(0) )
                chunk_size = int( len(res_list_unsliced) / n_chunk_left * 1.1)
                n_chunk_left -= 1

    return sliced_list_final
####################################
def pdb_slice_with_patching( input_pdb, out_name, slice_res_list ) :
    """
    Slice the pdb file with given input residue list. 
    Patch the nearby residues to the sliced pdb.
    Return a sorted new residue number list and a list of patched res (in new order).
    """
    check_path_exist( input_pdb ) 

    dist_cutoff = 5.0

    patched_res = find_nearby_res( input_pdb, slice_res_list, dist_cutoff )
    patched_res.sort()
    all_res = slice_res_list + patched_res
    all_res.sort()
    patched_res_new = []
    for res in patched_res :
        patched_res_new.append( all_res.index(res) + 1 )

    pdb_slice( input_pdb, out_name, all_res )
    return [all_res, patched_res_new]
####################################
def sliced2orig_merge_back( orig_pdb, new_pdb, out_name, res_list ) :
    """
    Merge processed sliced segment back to original pdb (Rosetta format).
    """
    check_path_exist( orig_pdb )
    check_path_exist( new_pdb )

    out = open(out_name, 'w')
    new_pdb_read = open(new_pdb)
    new_pdb_line = new_pdb_read.readline()
    res_new_pre = -10
    res_orig_pre = -10
    is_residue_done = False
    atom_num = 0
    for line in open(orig_pdb) :
        if len(line) < 4 :
            out.write(line)
            continue
        if line[0:4] == 'ATOM' :
            res_num = int(line[22:26])
            if res_num != res_orig_pre :
                res_orig_pre = res_num
                is_residue_done = False
            if is_residue_done :
                continue
            if not res_num in res_list :
                atom_num += 1
                out.write('%s%7d%s' %(line[0:4], atom_num, line[11:]) )
            else :
                if len(new_pdb_line) > 4 and new_pdb_line[0:4] == 'ATOM' :
                    atom_num += 1
                    out.write('%s%7d%s%4d%s' % (new_pdb_line[0:4], atom_num, new_pdb_line[11:22], res_num, new_pdb_line[26:]) )
                    res_new_pre = int( new_pdb_line[22:26] )

                while True :
                    new_pdb_line = new_pdb_read.readline()
                    if new_pdb_line == "" :
                        break
                    if len(new_pdb_line) > 4 and new_pdb_line[0:4] == 'ATOM' :
                        res_new = int( new_pdb_line[22:26] )
                        if res_new == res_new_pre :
                            atom_num += 1
                            out.write('%s%7d%s%4d%s' % (new_pdb_line[0:4], atom_num, new_pdb_line[11:22], res_num, new_pdb_line[26:]) )
                            res_new_pre = int( new_pdb_line[22:26] )
                        else :
                            break
                is_residue_done = True
    out.close()
    return True
