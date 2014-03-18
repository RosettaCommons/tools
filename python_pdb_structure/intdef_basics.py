
def write_subset( intres, chain ) :
    return write_subset_for_chain( intres[ chain ], chain )

def write_subset_for_chain( intres_for_chain, chain ) :
    lines = []
    for ires in intres_for_chain:
        lines.append( chain + " " + str( ires ) + "\n" )
    return lines

def load_subset( subset_file ):
    subset = []
    lines = open( subset_file ).readlines()
    for line in lines:
        if line == "\n" :
            continue
        subset.append( ( line.partition(" ")[ 0 ], line.partition(" ")[ 2 ].strip() ))
    return subset

def chains_from_xtal_chains_filename( pdbfile ) :
    pdb4char = pdbfile[ 0:4 ]
    chainA = pdbfile[ -6 ]
    chainB = pdbfile[ -5 ]
    return pdb4char, chainA, chainB

def intdef_fname_for_ssdes( pdb4char, chain ) :
    return pdb4char + "_ssintdef_" + chain + ".txt"

def intdef_fname_for_dsdes( pdb4char ):
    return pdb4char + "_dsintdef.txt"
