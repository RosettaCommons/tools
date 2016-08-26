#!/usr/bin/env python
import sys, os, os.path

def read_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    return my_file

def clear_old_files(fname):
    if os.path.isfile('%s.log' % fname) == True:
        os.remove('%s.log' % fname)
    if os.path.isfile('%s.cmdpath' % fname) == True:
        os.remove('%s.cmdpath' % fname)
    if os.path.isfile('%s.list' % fname) == True:
        os.remove('%s.list' % fname)

def write_path_files(block, fname):
    #do not clear, called from write_x_blocks
    with open('%s.cmdpath' % fname, 'a') as mypathfile:
        mypathfile.write('%s \n' % block[2])
    with open('%s.list' % fname, 'a') as mypathfile:
        mypathfile.write('%s \n' % block[0])

def write_path_files_from_blocks(block_set, fname):
    clear_old_files(fname)
    for block in block_set:
        write_path_files(block, fname)

def write_trim_blocks(block_set, fname):
    clear_old_files(fname)
    for block in block_set:
        new_block = []
        for line in block[1]:
            by_col = (''.join((line.lower()).split(':'))).split()
            if 'error' in by_col or '[error]' in by_col:
                new_block.append(line)
        with open('%s.log' % fname, 'a') as myfile:
            myfile.write('\n\n\n******%s******\n\n\n' % block[0])
            if len(block) == 4:
                myfile.write(block[3])
            myfile.writelines([('%s  {0}' % block[0]).format(i) for i in new_block])
        write_path_files(block, fname)

def write_full_blocks(block_set, fname):
    clear_old_files(fname)
    for block in block_set:
        with open('%s.log' % fname, 'a') as myfile:
            myfile.write('\n\n\n******%s******\n\n\n' % block[0])
            myfile.writelines([('%s  {0}' % block[0]).format(i) for i in block[1]])
        write_path_files(block, fname)

def log_blocker(log):
    log_blocks = []

    # points of interest (starts)
    indices = []
    for ind, line in enumerate(log):
        if line[:26] == "core.init: Rosetta version":
            indices.append(ind)
    # paired points based on real starts and assumed ends
    index_blocks = []
    for i, ind_num in enumerate(indices):
        if i == len(indices)-1:
            index_blocks.append([ind_num, len(log)])
        else:
            index_blocks.append([ind_num, indices[i+1]-1])
    print "The number of pdb log files via index_blocks: ", len(index_blocks), "\n"
    #print index_blocks
    count = 0
    # Log_blocks is a list of lists where each list has three members
    # member one is the pdb
    # member two is the block (chunk of log file)
    # member three is the path (exe + file path)
    log_blocks = []
    block = []
    for i, line in enumerate(log):
        if line[:20] == "core.init: command: ":
            path = line.split('command: ')[-1].strip('\r\n') + ' >>log 2>&1'
        if line[:51] == "protocols.jd2.PDBJobInputter: filling pose from PDB" or \
           line[:55] == "protocols.jd2.mmCIFJobInputter: filling pose from mmCIF":
            pdb = line.split('/')[-1][3:7]
        if i >= index_blocks[count][0] and i <= index_blocks[count][1]:
            block.append(line)

        if i == index_blocks[count][1]:
            log_blocks.append([pdb,block,path])
            block = []
            count += 1
        elif i == len(log)-1:
            log_blocks.append([pdb,block,path])
    print "The number of pdb log files via log_blocks dictated by index_blocks ", len(log_blocks), "\n"

    return log_blocks

def identify_blocks_by_pdb(block):
    for line in block:
        if line[:51] == "protocols.jd2.PDBJobInputter: filling pose from PDB":
            pdb = line.split('/')[-1][3:7]
            #print pdb
    return [pdb, block]

def identify_errors(log_blocks):
    # Each of these return data structures are like the log_blocks
    # however they have less information
    error_blocks = []                       #1
    fill_error_blocks = []                  #3
    ace_error_blocks = []                   #2
    assert_segfault_blocks = []             #17
    misc_segfault_blocks = []               #0
    resMap_range_error_blocks = []          #4
    partly_recognized_error_blocks = []     #5
    rotno_error_blocks = []                 #6
    polymer_bond_error_blocks = []          #7
    staple_error_blocks = []                #8
    aceCYS_error_blocks = []                #9
    unREC_res_error_blocks = []             #10
    unREC_ele_error_blocks = []             #11
    unREC_aType_error_blocks = []           #12
    unREC_token_error_blocks = []           #13
    unREC_expTec_error_blocks = []          #14
    pose_load_error_blocks = []             #15
    typeEle_error_blocks = []               #16
    zero_length_xyzVector_error_blocks = [] #18
    sugar_variant_error_blocks = []         #19
    zero_atom_restype_error_blocks = []     #20
    multiple_disulfides_error_blocks = []   #21
    seqpos_0_error_blocks = []              #22
    disulfide_restype_error_blocks = []     #23
    disulfide_from_atom_error_blocks = []   #24
    header_compound_value_error_blocks = [] #25
    multiple_nr_SS_error_blocks = []        #26
    zero_nres_blocks = []                   #27

    for block in log_blocks:
        pdb = block[0]
        path = block[2]
        block = block[1]
        assertion = False
        for ind, line in enumerate(block):
            #print line
            # splits line by colons and spaces
            by_col = (''.join((line.lower()).split(':'))).split()
            #by_braq = (line.lower()).split(']')
            if ind == len(block)-1:
                by_star = (line.lower()).split('*')
                if 'completed' not in by_star:
                    for line_num in xrange(len(block)-1):
                        rvs_line = (block[ind - line_num].lower()).split()
                        if "assertion" in rvs_line and 'failed.' in rvs_line: #[1].split():
                            assert_segfault_blocks.append([pdb,block,path,''.join(rvs_line)])
                            assertion = True
                            break
                    if assertion == False:
                        misc_segfault_blocks.append([pdb,block,path])
            if 'error' in by_col:
                # For each error line we break up the +- 1 lines by spaces
                next_line = (block[ind+1].lower()).split()
                pre_line = (block[ind-1].lower()).split()
                if 'ace' in next_line:
                    ace_error_blocks.append([pdb,block,path])
                    break
                elif "fill_missing_atoms!" in by_col:#[1].split():
                    fill_error_blocks.append([pdb,block,path])
                    break
                elif "packed_rotno_conversion_data_current_" in by_col:#[1].split():
                    rotno_error_blocks.append([pdb,block,path])
                    break
                elif "polymer" in by_col and 'incompatible' in by_col: #[1].split():
                    polymer_bond_error_blocks.append([pdb,block,path])
                    break
                elif "patchoperation" in by_col: #[1].split():
                    staple_error_blocks.append([pdb,block,path])
                    break
                elif "disulfide-bonded" in by_col: #and "partner" in by_col[1].split():
                    aceCYS_error_blocks.append([pdb,block,path])
                    break
                elif '3-letter' in next_line:
                    partly_recognized_error_blocks.append([pdb,block,path,' '.join(next_line) + '\n'])
                    break
                elif "unrecognized" in by_col and "residue" in by_col: #[1].split():
                    unREC_res_error_blocks.append([pdb,block,path])
                    break
                elif "unrecognized" in by_col and "element_symbol" in by_col: #[1].split():
                    unREC_ele_error_blocks.append([pdb,block,path])
                    break
                elif "unrecognized" in by_col and "atom_type_name" in by_col: #[1].split():
                    unREC_aType_error_blocks.append([pdb,block,path])
                    break
                elif "unrecognized" in pre_line and "compound" in pre_line and "token" in pre_line: #[1].split():
                    unREC_token_error_blocks.append([pdb,block,path])
                    break
                elif "unrecognized" in pre_line and "experimental" in pre_line and "technique" in pre_line: #[1].split():
                    unREC_expTec_error_blocks.append([pdb,block,path])
                    break
                elif "attempting" in pre_line and "compound" in pre_line and "value" in pre_line: #[1].split():
                    header_compound_value_error_blocks.append([pdb,block,path])
                    break
                elif 'res_map' in next_line and 'range' in next_line:
                    resMap_range_error_blocks.append([pdb,block,path])
                    break
                elif "normalize" in next_line and "xyzvector" in next_line and "zero" in next_line:
                    zero_length_xyzVector_error_blocks.append([pdb,block,path])
                    break
                #unclear that pose_load_error is real - I think this catches when a file is missing?
                elif "exception" in by_col and "jobdistributor" in by_col: #[1].split():
                    pose_load_error_blocks.append([pdb,block,path])
                    break
                elif "cannot" in by_col and "type" in by_col and "element" in by_col: #[1].split():
                    typeEle_error_blocks.append([pdb,block,path])
                    break
                elif "unable to find desired variant residue" in line:
                    sugar_variant_error_blocks.append([pdb,block,path])
                    break
                elif "Cannot load in ResidueType for entry with no atoms." in line:
                    zero_atom_restype_error_blocks.append([pdb,block,path])
                    break
                elif "SSBond records list multiple disulfides for this residue!" in line:
                    multiple_disulfides_error_blocks.append([pdb,block,path])
                    break
                elif "SSBond records list multiple nonredundant disulfides for this residue!" in line:
                    multiple_nr_SS_error_blocks.append([pdb,block,path])
                    break
                elif "The sequence position requested was 0." in line:
                    seqpos_0_error_blocks.append([pdb,block,path])
                    break
                elif "unable to create appropriate residue type for disulfide" in line:
                    disulfide_restype_error_blocks.append([pdb,block,path])
                    break
                elif "Can't find an atom to disulfide bond from" in line:
                    disulfide_from_atom_error_blocks.append([pdb,block,path])
                    break
                else:
                    error_blocks.append([pdb,block,path])
                    break
            elif '[error]' in by_col:
                if "Cannot normalize xyzVector of length() zero" in line: #this one has a race condition sometimes, shows up on the wrong line
                    zero_length_xyzVector_error_blocks.append([pdb,block,path])
                    break
                else:
                    error_blocks.append([pdb,block,path])
                    break
            elif "nres 0" in line:
                zero_nres_blocks.append([pdb,block,path])
                break

    return [misc_segfault_blocks, \
            error_blocks, \
            ace_error_blocks, \
            fill_error_blocks, \
            resMap_range_error_blocks, \
            partly_recognized_error_blocks, \
            rotno_error_blocks, \
            polymer_bond_error_blocks, \
            staple_error_blocks, \
            aceCYS_error_blocks, \
            unREC_res_error_blocks   , \
            unREC_ele_error_blocks   , \
            unREC_aType_error_blocks , \
            unREC_token_error_blocks , \
            unREC_expTec_error_blocks, \
            pose_load_error_blocks, \
            typeEle_error_blocks, \
            assert_segfault_blocks, \
            zero_length_xyzVector_error_blocks, \
            sugar_variant_error_blocks, \
            zero_atom_restype_error_blocks, \
            multiple_disulfides_error_blocks, \
            seqpos_0_error_blocks, \
            disulfide_restype_error_blocks, \
            disulfide_from_atom_error_blocks, \
            header_compound_value_error_blocks, \
            multiple_nr_SS_error_blocks, \
            zero_nres_blocks] #27
def main(argv):

    log_file = read_file(argv[0])

    log_blocks = log_blocker(log_file)

    #log_blocks = [identify_blocks_by_pdb(x) for x in log_blocks]

    all_errors = identify_errors(log_blocks)

    print len(all_errors[ 0]), "structures with misc segfaults [often misclassified due to output buffering between cerr and cout] (miscSegfault.log) " , len(all_errors[ 0])
    print len(all_errors[17]), "structures with null pointer assertions (nullPointerAssertion.log) ", len(all_errors[17])
    print len(all_errors[ 1]), "structures with misc unidentified errors (unidentified_error) ", len(all_errors[1])
    print len(all_errors[ 2]), "structures with acetylated N-terminus errors (ACE_error.log) ", len(all_errors[2])
    print len(all_errors[ 3]), "structures with 'too many tries in fill_missing_atoms!' errors (fill_missing_atoms.log) ", len(all_errors[3])
    print len(all_errors[ 4]), "structures with 'Residue outside res_map range' errors (resMap_range_error.log) ", len(all_errors[4])
    print len(all_errors[ 5]), "structures with partly recognized residue errors (usually branch-point) (partly_recognized_res.log)", len(all_errors[5])
    print len(all_errors[ 6]), "structures with packed_rotno_conversion_data_current_ errors (rotno_error.log) ", len(all_errors[6])
    print len(all_errors[ 7]), "structures with 'Can't create a polymer bond' errors (usually sugars) (polymer_bond_error.log) ", len(all_errors[7])
    print len(all_errors[ 8]), "structures with PatchOperation errors (PatchOperation.log) ", len(all_errors[8])
    print len(all_errors[ 9]), "structures with ace.CYS errors (needs acetylated CYS restype) (ace_CYS.log) ", len(all_errors[9])
    print len(all_errors[10]), "structures with unrecognized residue errors (straight 'unrecognized residue', usually ligands?) (resUnrec.log)  ", len(all_errors[10])
    print len(all_errors[11]), "structures with unrecognized element errors (eleUnrec.log) ", len(all_errors[11])
    print len(all_errors[12]), "structures with unrecognized atom_type_name errors (aTypeUnrec.log) ", len(all_errors[12])
    print "next few are header errors"
    print len(all_errors[13]), "structures with 'unrecognized compound token string' errors (token.log) ", len(all_errors[13])
    print len(all_errors[14]), "structures with unrecognized experimental technique errors (expTech.log) ", len(all_errors[14])
    print len(all_errors[25]), "structures with compound header errors (NOTE TYPO IN ROSETTA LOG!) (compoundHeader.log) ", len(all_errors[25])
    print len(all_errors[15]), "structures with misc pose load errors (usually means missing file on my end) (poseLoad.log) ", len(all_errors[15])
    print len(all_errors[16]), "structures with cannot type atom with element errors (typAtomEle.log) ", len(all_errors[16])
    print len(all_errors[18]), "structures with 'Cannot normalize xyzVector of length() zero' errors (zeroLengthXYZVector.log) ", len(all_errors[18])
    print len(all_errors[19]), "structures with 'unable to find desired variant residue' errors [usually sugars]' (sugarVariant.log) ", len(all_errors[19])
    print len(all_errors[20]), "structures with 'Cannot load in ResidueType for entry with no atoms.' (0AtomRestype.log) ", len(all_errors[20])
    print len(all_errors[21]), "structures with 'SSBond records list multiple disulfides for this residue.' (multiSS.log) ", len(all_errors[21])
    print len(all_errors[26]), "structures with 'SSBond records list multiple nonredundant disulfides for this residue.' (multiNR_SS.log) ", len(all_errors[26])
    print len(all_errors[22]), "structures with 'The sequence position requested was 0.' (seqpos0.log) ", len(all_errors[22])
    print len(all_errors[23]), "structures with 'unable to create appropriate residue type for disulfide' (resTypeSS.log) ", len(all_errors[23])
    print len(all_errors[24]), "structures with 'Can't find an atom to disulfide bond from' (atomSS.log) ", len(all_errors[24])
    print len(all_errors[27]), "structures with zero residue structures' (zeronres.list) ", len(all_errors[27])

    write_full_blocks(all_errors[ 0], 'miscSegfault')
    write_full_blocks(all_errors[ 1], 'unidentified_error')
    write_trim_blocks(all_errors[ 2], 'ACE_error')
    write_trim_blocks(all_errors[ 3], 'fill_missing_atoms')
    write_trim_blocks(all_errors[ 4], 'resMap_range_error')
    write_trim_blocks(all_errors[ 5], 'partly_recognized_res')
    write_trim_blocks(all_errors[ 6], 'rotno_error')
    write_trim_blocks(all_errors[ 7], 'polymer_bond_error')
    write_trim_blocks(all_errors[ 8], 'PatchOperation')
    write_trim_blocks(all_errors[ 9], 'ace_CYS')
    write_trim_blocks(all_errors[10], 'resUnrec')
    write_trim_blocks(all_errors[11], 'eleUnrec')
    write_trim_blocks(all_errors[12], 'aTypeUnrec')
    write_full_blocks(all_errors[13], 'token')
    write_full_blocks(all_errors[14], 'expTech')
    write_trim_blocks(all_errors[15], 'poseLoad')
    write_trim_blocks(all_errors[16], 'typAtomEle')
    write_full_blocks(all_errors[17], 'nullPointerAssertion')
    write_trim_blocks(all_errors[18], 'zeroLengthXYZVector')
    write_trim_blocks(all_errors[19], 'sugarVariant')
    write_trim_blocks(all_errors[20], '0AtomRestype')
    write_trim_blocks(all_errors[21], 'multiSS')
    write_trim_blocks(all_errors[22], 'seqpos0')
    write_trim_blocks(all_errors[23], 'resTypeSS')
    write_trim_blocks(all_errors[24], 'atomSS')
    write_trim_blocks(all_errors[25], 'compoundHeader')
    write_trim_blocks(all_errors[26], 'multiNR_SS')
    write_path_files_from_blocks(all_errors[27], 'zeronres')

if __name__ == '__main__':

    main(sys.argv[1:])
