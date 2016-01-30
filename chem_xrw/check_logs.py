#!/usr/bin/python
import sys, os, os.path

def read_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    return my_file

def write_trim_blocks(block_set, fname):
    if os.path.isfile(fname) == True:
        os.remove(fname)
    for block in block_set:
        new_block = []
        for line in block[1]:
            by_col = (''.join((line.lower()).split(':'))).split()
            if 'error' in by_col or '[error]' in by_col:
                new_block.append(line)
        with open(fname, 'a') as myfile:
            myfile.write('\n\n\n******%s******\n\n\n' % block[0])
            if len(block_set) == 3:
                myfile.write(block[2])
            myfile.writelines([('%s  {0}' % block[0]).format(i) for i in new_block])
            
def write_full_blocks(block_set, fname):
    if os.path.isfile(fname) == True:
        os.remove(fname)
    for block in block_set:
        with open(fname, 'a') as myfile:
            myfile.write('\n\n\n******%s******\n\n\n' % block[0])
            myfile.writelines([('%s  {0}' % block[0]).format(i) for i in block[1]])

def log_blocker(log):
    log_blocks = []
    
    indices = []
    for ind, line in enumerate(log):
        if line[:26] == "core.init: Rosetta version":
            indices.append(ind)
    #print len(indices)
    index_blocks = []
    for i, ind_num in enumerate(indices):
        if i == len(indices)-1:
            index_blocks.append([ind_num, len(log)])
        else:
            index_blocks.append([ind_num, indices[i+1]-1])
    print "The number of pdb log files via index_blocks: ", len(index_blocks), "\n"
    #print index_blocks
    count = 0
    log_blocks = []
    block = []
    for i, line in enumerate(log):
        if line[:51] == "protocols.jd2.PDBJobInputter: filling pose from PDB":
            pdb = line.split('/')[-1][3:7]
        if i >= index_blocks[count][0] and i <= index_blocks[count][1]:
            block.append(line)
            
        if i == index_blocks[count][1]:
            log_blocks.append([pdb,block])
            block = []
            count += 1
        elif i == len(log)-1:
            log_blocks.append([pdb,block])
    print "The number of pdb log files via log_blocks dictated by index_blocks ", len(log_blocks), "\n"
    #print [x[0] for x in log_blocks]
    
    #for iblock in index_blocks:
        #log_block = [line for i, line in enumerate(log) if i >= iblock[0] and i <= iblock[1]]
        #log_blocks.append(log_block)
    
    #print log_blocks[0]
    return log_blocks, len(index_blocks) 
            
def identify_blocks_by_pdb(block):
    for line in block:
        if line[:51] == "protocols.jd2.PDBJobInputter: filling pose from PDB":
            pdb = line.split('/')[-1][3:7]
            #print pdb
    return [pdb, block]

def identify_errors(log_blocks):
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

    for block in log_blocks:
        pdb = block[0]
        block = block[1]
        assertion = False
        for ind, line in enumerate(block):
            #print line
            #by_col = (line.lower()).split(':')
            by_col = (''.join((line.lower()).split(':'))).split()
            #by_braq = (line.lower()).split(']')
            if ind == len(block)-1:
                by_star = (line.lower()).split('*')
                if 'completed' not in by_star:
                    for line_num in xrange(len(block)-1):
                        rvs_line = (block[ind - line_num].lower()).split()
                        if "assertion" in rvs_line and 'failed.' in rvs_line: #[1].split():
                            assert_segfault_blocks.append([pdb,block,''.join(rvs_line)])
                            assertion = True
                            break
                    if assertion == False:
                        misc_segfault_blocks.append([pdb,block])
            if 'error' in by_col:
                next_line = (block[ind+1].lower()).split()
                pre_line = (block[ind-1].lower()).split()
                if 'ace' in next_line:
                    ace_error_blocks.append([pdb,block])
                    break
                elif "fill_missing_atoms!" in by_col:#[1].split():
                    fill_error_blocks.append([pdb,block])
                    break
                elif "packed_rotno_conversion_data_current_" in by_col:#[1].split():
                    rotno_error_blocks.append([pdb,block])
                    break
                elif "polymer" in by_col and 'incompatible' in by_col: #[1].split():
                    polymer_bond_error_blocks.append([pdb,block])
                    break
                elif "patchoperation" in by_col: #[1].split():
                    staple_error_blocks.append([pdb,block])
                    break
                elif "disulfide-bonded" in by_col: #and "partner" in by_col[1].split():
                    aceCYS_error_blocks.append([pdb,block])
                    break
                elif '3-letter' in next_line:
                    letter3_error_blocks.append([pdb,block,''.join(next_line)])
                    break
                elif "unrecognized" in by_col and "residue" in by_col: #[1].split():
                    unREC_res_error_blocks.append([pdb,block])
                    break
                elif "unrecognized" in by_col and "element_symbol" in by_col: #[1].split():
                    unREC_ele_error_blocks.append([pdb,block])
                    break
                elif "unrecognized" in by_col and "atom_type_name" in by_col: #[1].split():
                    unREC_aType_error_blocks.append([pdb,block])
                    break
                elif "unrecognized" in pre_line and "compound" in pre_line and "token" in pre_line: #[1].split():
                    unREC_token_error_blocks.append([pdb,block])
                    break
                elif "unrecognized" in pre_line and "experimental" in pre_line and "technique" in pre_line: #[1].split():
                    unREC_expTec_error_blocks.append([pdb,block])
                    break
                elif 'res_map' in next_line and 'range' in next_line:
                    resMap_range_error_blocks.append([pdb,block])
                    break
                elif "exception" in by_col and "jobdistributor" in by_col: #[1].split():
                    pose_load_error_blocks.append([pdb,block])
                    break
                elif "cannot" in by_col and "type" in by_col and "element" in by_col: #[1].split():
                    typeEle_error_blocks.append([pdb,block])
                    break
                else:
                    error_blocks.append([pdb,block])
                    break
            elif '[error]' in by_col:
                error_blocks.append([pdb,block])
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
            assert_segfault_blocks]


def main(argv):

    log_file = read_file(argv[0])
    
    log_blocks, pdb_per_block = log_blocker(log_file)
    
    #log_blocks = [identify_blocks_by_pdb(x) for x in log_blocks]

    all_errors = identify_errors(log_blocks)

    print "The number of blocks with misc segfaults ", len(all_errors[0]), "\n"
    
    print "The number of blocks with assertion segfaults ", len(all_errors[17]), "\n"
    
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
    
    print "The number of blocks with pose load errors ", len(all_errors[15]), "\n"

    print "The number of blocks with cannot type atom with element errors ", len(all_errors[16]), "\n"
    
    write_full_blocks(all_errors[ 0], 'miscSegfault.log')
    write_trim_blocks(all_errors[ 1], 'unidentified_error.log')
    write_trim_blocks(all_errors[ 2], 'ACE_error.log')
    write_trim_blocks(all_errors[ 3], 'fill_atom_error.log')
    write_trim_blocks(all_errors[ 4], 'resMap_range_error.log')
    write_trim_blocks(all_errors[ 5], 'nonACE_res_error.log')
    write_trim_blocks(all_errors[ 6], 'rotno_error.log')
    write_trim_blocks(all_errors[ 7], 'polymer_bond_error.log')
    write_trim_blocks(all_errors[ 8], 'PatchOperation_error.log')
    write_trim_blocks(all_errors[ 9], 'ace_CYS.log')
    write_trim_blocks(all_errors[10], 'resUnrec.log')
    write_trim_blocks(all_errors[11], 'eleUnrec.log')
    write_trim_blocks(all_errors[12], 'aTypeUnrec.log')
    write_trim_blocks(all_errors[13], 'token.log')
    write_trim_blocks(all_errors[14], 'expTech.log')
    write_trim_blocks(all_errors[15], 'poseLoad.log')
    write_trim_blocks(all_errors[16], 'typAtomEle.log')
    write_trim_blocks(all_errors[17], 'assertSegfault.log')


if __name__ == '__main__':

    main(sys.argv[1:])
