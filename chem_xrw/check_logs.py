import sys, os

def read_file(file_path):
    with open(file_path, 'r') as f:
        my_file = f.readlines()
    return my_file


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
    print len(index_blocks)
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
    print len(log_blocks)
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
    segfault_blocks = []
    for block in log_blocks:
        pdb = block[0]
        block = block[1]
        for ind, line in enumerate(block):
            #print line
            by_col = (line.lower()).split(':')
            by_braq = (line.lower()).split(']')
            if ind == len(block)-1:
                by_star = (line.lower()).split('*')
                if 'completed' not in by_star:
                    segfault_blocks.append([pdb,block])
            if 'error' in by_col:
                error_blocks.append([pdb,block])
                break
            elif '[error' in by_braq:
                error_blocks.append([pdb,block])
                break
    return error_blocks, segfault_blocks

def main(argv):

    log_file = read_file(argv[0])
    
    log_blocks, pdb_per_block = log_blocker(log_file)
    
    #log_blocks = [identify_blocks_by_pdb(x) for x in log_blocks]

    error_blocks, segfault_blocks = identify_errors(log_blocks)

    print len(error_blocks)
    
    print len(segfault_blocks)
        

if __name__ == '__main__':

    main(sys.argv[1:])
