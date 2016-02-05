#!/usr/bin/python
import sys, os, os.path

def read_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    return my_file

def write_trim_blocks(block_set, fname):
    if os.path.isfile('%s.log' % fname) == True:
        os.remove('%s.log' % fname)
    if os.path.isfile('%s.path' % fname) == True:
        os.remove('%s.path' % fname)
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
        with open('%s.path' % fname, 'a') as mypathfile:
            mypathfile.write('%s \n' % block[2])
            
def write_full_blocks(block_set, fname):
    if os.path.isfile('%s.log' % fname) == True:
        os.remove('%s.log' % fname)
    if os.path.isfile('%s.path' % fname) == True:
        os.remove('%s.path' % fname)
    for block in block_set:
        with open('%s.log' % fname, 'a') as myfile:
            myfile.write('\n\n\n******%s******\n\n\n' % block[0])
            myfile.writelines([('%s  {0}' % block[0]).format(i) for i in block[1]])
        with open('%s.path' % fname, 'a') as mypathfile:
            mypathfile.write('%s \n' % block[2])

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
    error_blocks = []
    fill_error_blocks = []
    ace_error_blocks = []
    assert_segfault_blocks = []
    misc_segfault_blocks = []
    resMap_range_error_blocks = []
    letter3_error_blocks = []
    rotno_error_blocks = []
    polymer_bond_error_blocks = []
    staple_error_blocks = []
    aceCYS_error_blocks = []
    unREC_res_error_blocks = [] 
    unREC_ele_error_blocks = [] 
    unREC_aType_error_blocks = [] 
    unREC_token_error_blocks = [] 
    unREC_expTec_error_blocks = [] 
    pose_load_error_blocks = []
    typeEle_error_blocks = []
    zero_length_xyzVector_error_blocks = []

    for block in log_blocks:
        pdb = block[0]
        path = block[2]
        block = block[1]
        assertion = False
        for ind, line in enumerate(block):
            #print line
            # splits line by a colon 
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
                    letter3_error_blocks.append([pdb,block,path,' '.join(next_line) + '\n'])
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
                elif 'res_map' in next_line and 'range' in next_line:
                    resMap_range_error_blocks.append([pdb,block,path])
                    break
                elif "cannot normalize xyzvector of length() zero" in by_col:
                    zero_length_xyzVector_error_blocks.append([pdb,block,path])
                    break
                    #unclear that pose_load_error is real - I think this catches when a file is missing?
                elif "exception" in by_col and "jobdistributor" in by_col: #[1].split():
                    pose_load_error_blocks.append([pdb,block,path])
                    break
                elif "cannot" in by_col and "type" in by_col and "element" in by_col: #[1].split():
                    typeEle_error_blocks.append([pdb,block,path])
                    break
                else:
                    error_blocks.append([pdb,block,path])
                    break
            elif '[error]' in by_col:
                error_blocks.append([pdb,block,path])
                break

    return [misc_segfault_blocks, \
            error_blocks, \
            ace_error_blocks, \
            fill_error_blocks, \
            resMap_range_error_blocks, \
            letter3_error_blocks, \
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
            zero_length_xyzVector_error_blocks]


def main(argv):

    log_file = read_file(argv[0])
    
    log_blocks = log_blocker(log_file)
    
    #log_blocks = [identify_blocks_by_pdb(x) for x in log_blocks]

    all_errors = identify_errors(log_blocks)

    print "The number of blocks with misc segfaults ", len(all_errors[0]), "\n"
    
    print "The number of blocks with null pointer assertions ", len(all_errors[17]), "\n"
    
    print "The number of blocks with misc unidentified errors ", len(all_errors[1]), "\n"

    print "The number of blocks with ace errors ", len(all_errors[2]), "\n"

    print "The number of blocks with fill errors ", len(all_errors[3]), "\n"
    
    print "The number of blocks with resMap range errors ", len(all_errors[4]), "\n"
    
    print "The number of blocks with res errors not ace ", len(all_errors[5]), "\n"
        
    print "The number of blocks with rotno errors ", len(all_errors[6]), "\n"
    
    print "The number of blocks with polymer bond errors ", len(all_errors[7]), "\n"
    
    print "The number of blocks with PatchOperation errors ", len(all_errors[8]), "\n"
    
    print "The number of blocks with ace.CYS errors ", len(all_errors[9]), "\n"
    
    print "The number of blocks with unrecognized residue errors ", len(all_errors[10]), "\n"
    
    print "The number of blocks with unrecognized element errors ", len(all_errors[11]), "\n"
    
    print "The number of blocks with unrecognized atom_type errors ", len(all_errors[12]), "\n"
    
    print "The number of blocks with unrecognized token errors ", len(all_errors[13]), "\n"
    
    print "The number of blocks with unrecognized experimental_technique errors ", len(all_errors[14]), "\n"
    
    print "The number of blocks with 'Cannot normalize xyzVector of length() zero' errors", len(all_errors[17])

    print "The number of blocks with pose load errors ", len(all_errors[15]), "\n"

    print "The number of blocks with cannot type atom with element errors ", len(all_errors[16]), "\n"
    
    print "The number of blocks with 'Cannot normalize xyzVector of length() zero' errors", len(all_errors[17])

    write_full_blocks(all_errors[ 0], 'miscSegfault')
    write_trim_blocks(all_errors[ 1], 'unidentified_error')
    write_trim_blocks(all_errors[ 2], 'ACE_error')
    write_trim_blocks(all_errors[ 3], 'fill_atom_error')
    write_trim_blocks(all_errors[ 4], 'resMap_range_error')
    write_trim_blocks(all_errors[ 5], 'nonACE_res_error')
    write_trim_blocks(all_errors[ 6], 'rotno_error')
    write_trim_blocks(all_errors[ 7], 'polymer_bond_error')
    write_trim_blocks(all_errors[ 8], 'PatchOperation_error')
    write_trim_blocks(all_errors[ 9], 'ace_CYS')
    write_trim_blocks(all_errors[10], 'resUnrec')
    write_trim_blocks(all_errors[11], 'eleUnrec')
    write_trim_blocks(all_errors[12], 'aTypeUnrec')
    write_trim_blocks(all_errors[13], 'token')
    write_trim_blocks(all_errors[14], 'expTech')
    write_trim_blocks(all_errors[15], 'poseLoad')
    write_trim_blocks(all_errors[16], 'typAtomEle')
    write_full_blocks(all_errors[17], 'nullPointerAssertion')


if __name__ == '__main__':

    main(sys.argv[1:])
