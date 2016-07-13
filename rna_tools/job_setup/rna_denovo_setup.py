#!/usr/bin/env python

from rna_server_conversions import prepare_fasta_and_params_file_from_sequence_and_secstruct, get_all_stems
from sys import argv
from os import system, getcwd
from os.path import exists, dirname, basename
from parse_options import parse_options, get_ints, get_resnum_chain
from make_tag import make_tag, make_tag_with_dashes
import string
from rosetta_exe import rosetta_exe

# I could make this a little bit smarter:
# if I look at input_res, I should be able to figure out obligate pairs "on the fly", right? This
#  could either happen here, or within Rosetta.

def flatten(l):
    new_l = []
    for e in l:
        new_l.extend(e)
    return new_l

def Help():
    print
    print argv[0] + " 'sequence in quotes'  'secstruct in quotes'  [tag]"
    print
    print "example:  "+argv[0]+" 'aaa uuu'  '((( )))' test"
    print

if len(argv) < 3:
    Help()
    exit()


argv_input = []
for arg in argv: argv_input.append( arg ) # deepcopy

# note -- following has not been converted over entirely to handle input of residues through residue-numbering & chains
out_script = parse_options( argv, "out_script", "README_FARFAR" )
nstruct = parse_options(argv, "nstruct", 500)
sequence = parse_options( argv, "sequence", "" )
secstruct = parse_options( argv, "secstruct", "")
secstruct_general = parse_options( argv, "secstruct_general", "")
(working_res,working_chain) = parse_options( argv, "working_res", [[-1],[" "]] )
offset = parse_options( argv, "offset", 0 )
tag = parse_options( argv, 'tag', "" )
fasta = parse_options( argv, 'fasta', "" )
secstruct_file = parse_options( argv, 'secstruct_file', "" )
input_pdbs = parse_options( argv, 's', [""] )
input_silent_files = parse_options( argv, 'silent', [""] )
(input_silent_res,input_silent_chain) = parse_options( argv, 'input_silent_res', [[-1],[" "]] )
fixed_stems = parse_options( argv, 'fixed_stems', False )
no_minimize = parse_options( argv, 'no_minimize', False )
is_cst_gap = parse_options( argv, 'cst_gap', False )
native_pdb = parse_options( argv, 'native', "" )
working_native_pdb = parse_options( argv, 'working_native', "" )
cst_file = parse_options( argv, 'cst_file', "" )
data_file = parse_options( argv, 'data_file', "" )
cutpoint_open_res_chain = parse_options( argv, "cutpoint_open", [[-1],[" "]] )
cutpoint_closed_res_chain = parse_options( argv, "cutpoint_closed", [[-1],[" "]] )
extra_minimize_res_chain = parse_options( argv, "extra_minimize_res", [[-1],[" "]] )
extra_minimize_chi_res_chain = parse_options( argv, "extra_minimize_chi_res",[[-1],[" "]]  )
virtual_anchor_res_chain = parse_options( argv, "virtual_anchor", [[-1],[" "]] )
obligate_pair_res_chain = parse_options( argv, "obligate_pair", [[-1],[" "]] )
obligate_pair = zip( obligate_pair_res_chain[0],obligate_pair_res_chain[1] )
obligate_pair_explicit = parse_options( argv, "obligate_pair_explicit", [""] )
remove_obligate_pair = parse_options( argv, "remove_obligate_pair", [-1] )
remove_pair = parse_options( argv, "remove_pair", [-1] )
chain_connection = parse_options( argv, "chain_connection", [""] )
bps_moves = parse_options( argv, "bps_moves", False )
rosetta_folder = parse_options( argv, "rosetta_folder", "" )
extension = parse_options( argv, "extension", "" )

system( 'rm -rf %s' % out_script ) # output file with Rosetta command line -- will be replaced by this script.
#input_res and cutpoint_closed changes to be auto-generated

#print argv
#assert( len( argv ) == 1 )

extra_args = ""
if len( argv ) > 1: extra_args = string.join( argv[1:] )

if len(secstruct_general) > 0 and not bps_moves:
    raise ValueError("cannot supply secstruct_general without bps_moves")

def is_even( num ):
    return (  2 * (num/2) == num ) # even

def working_res_map( vector, working_res, working_chain = [] ):
    if len( working_res ) == 0: return vector
    if len( vector ) == 0: return vector
    working_vector = []
    if isinstance( vector[ 0 ], list ): # its a list of residue, then list of chain
        assert( isinstance( vector[ 1 ], list ) )
        assert( len( working_chain ) == len( working_res ) )
        working_res_chain = zip( working_res, working_chain )
        for m in zip( vector[0], vector[1] ):
            if m in working_res_chain: working_vector.append( working_res_chain.index( m )+1 )
    else:
        for m in vector:
            if m in working_res: working_vector.append( working_res.index( m )+1 )

    return working_vector

full_model_res    = []
full_model_chain = []
if len( fasta ) > 0 :
    assert( len(sequence) == 0 )
    lines = open( fasta ).readlines()
    for line in lines: # worst fasta reader in the world.
        line = line.strip()
        if len( line ) > 0 and line[0] == '>':
            if len(sequence) > 0:
                cutpoint_open_res_chain[ 0 ].append( full_model_res[ -1 ] )
                cutpoint_open_res_chain[ 1 ].append( full_model_chain[ -1 ] )
            for fasta_tag in line.split()[1:]:
                get_resnum_chain( fasta_tag, full_model_res, full_model_chain )
            continue
        sequence += line

if len( full_model_res ) > 0 and len( full_model_res ) < len( sequence )+5:
    print 'Thought we found a tag supplying model res and chain, but not enough to match sequence ... assuming residue numbering is 1,2,...'
    full_model_res = []
if len( full_model_res ) == 0: full_model_res = range( offset+1, offset+len(sequence)+1)
assert( len( full_model_res ) == len( sequence ) )

if len( secstruct_file ) > 0 :
    assert( len(secstruct) == 0 )
    secstruct = open( secstruct_file ).readlines()[0].strip().replace( '\n','').replace('\r','')

if len( secstruct ) == 0:
    for m in range(len(sequence)): secstruct += '.'
if len( secstruct_general ) == 0:
    for m in range(len(sequence)): secstruct_general += '.'

if( not len( sequence ) == len( secstruct )):
    print sequence
    print secstruct
    print 'Length of sequence & secstruct do not match: ', len( sequence ), len( secstruct )
    exit( 1 )
assert( secstruct.count('(') == secstruct.count(')') )
assert( secstruct.count('[') == secstruct.count(']') )
assert( secstruct.count('{') == secstruct.count('}') )
if ( tag == '' ):
    tag = basename(getcwd())

assert( is_even( len(remove_pair) ) )
for m in remove_pair:
    pos = full_model_res.index( m )
    assert( pos > -1 and pos < len( secstruct) and secstruct[ pos ] != '.' )
    secstruct = secstruct[:pos] + '.' + secstruct[ pos+1:]

if len(working_res) == 0:
    working_res = full_model_res
    working_chain = full_model_chain

working_sequence = ''
working_secstruct = ''
working_secstruct_general = ''
#working_res.sort()
for i in range(len(working_res)):
    m = working_res[ i ]
    if i > 0 and ( ( m > working_res[ i-1 ] + 1 ) or ( len( working_chain ) > 0 and working_chain [ i ] != working_chain[ i-1 ] ) ):
        working_sequence  += ' '
        working_secstruct += ' '
        working_secstruct_general += ' '
    working_sequence  += sequence[  full_model_res.index( m ) ]
    working_secstruct += secstruct[ full_model_res.index( m ) ]
    working_secstruct_general += secstruct_general[ full_model_res.index( m ) ]

print working_sequence
print working_secstruct
print working_secstruct_general

assert( working_secstruct.count('(') == working_secstruct.count(')') )
assert( working_secstruct.count('[') == working_secstruct.count(']') )
####################################################
#Add by FCC: process PDBs and generate reasonable obligate pairs
def get_seq_and_resnum( pdb ):
    resnum = []
    chain = []
    seq = ''
    old_res = None
    old_chain = None
    for line in open(pdb):
        if len(line) > 50 and line[:4] == 'ATOM':
            curr_res = int( line[22:26] )
            curr_chain = line[21]
            curr_seq = line[19]
            if curr_res != old_res or curr_chain != old_chain:
                if ( (curr_res,curr_chain) in zip( resnum,chain ) ):
                        raise ValueError('Residue/chain numbers in %s is not unique!!' % pdb)
                resnum.append(curr_res)
                chain.append(curr_chain)
                old_res = curr_res
                old_chain = curr_chain
                seq += curr_seq
    return seq, resnum, chain

def get_silent_seq( silent ):
    for line in open(silent):
        if 'SEQUENCE' in line:
            return line.split()[1].upper()
    raise ValueError('No sequence info found in the input silent file!')

input_res = []
input_chain = []
resnum_list = []
chain_list = []
for pdb in input_pdbs:
    pdb_seq, resnum, chain = get_seq_and_resnum(pdb)
    pdb_seq = pdb_seq.lower()
    actual_seq = ''
    for i in resnum:
        if i in input_res: print('WARNING! Input residue %s exists in two pdb files!!' % i)
        actual_seq += sequence[ full_model_res.index( i )]
    if pdb_seq != actual_seq.lower():
        print pdb_seq
        print actual_seq
        raise ValueError('The sequence in %s does not match input sequence!!' % pdb)
    resnum_list.append(resnum)
    chain_list.append(chain)
    input_res   += resnum
    input_chain += chain

for silent in input_silent_files:
    actual_seq = ''
    seq = get_silent_seq( silent )
    len_seq = len(seq)
    if len_seq > len(input_silent_res):
        raise ValueError('input_silent_res is shorter then the total length of silent files!')
    resnum = input_silent_res[:len_seq]
    input_silent_res = input_silent_res[len_seq:]
    chain = input_silent_chain[:len_seq]
    input_silent_chain = input_silent_chain[len_seq:]
    for i in resnum:
        if i in input_res: print('WARNING! Input residue %s exists in two pdb files!!' % i)
        actual_seq += sequence[full_model_res.index(i)]
    if seq.lower() != actual_seq.lower():
        raise ValueError('The sequence in %s does not match input sequence!!' % silent)
    resnum_list.append(resnum)
    chain_list.append(chain)
    input_res += resnum
    input_chain += chain

def already_listed_in_obligate_pair( new_pair, obligate_pair ):
    new_pos1 = new_pair[ 0 ]
    new_pos2 = new_pair[ 1 ]
    already_listed = False
    for m in range( len( obligate_pair)/2 ):
        pos1 = obligate_pair[ 2*m ]
        pos2 = obligate_pair[ 2*m+1 ]
        if ( pos1 == new_pos1 and pos2 == new_pos2 ):
            already_listed = True
            break
    return already_listed

for resnum,chain in zip(resnum_list,chain_list):
    #Find obligate pairs
    chunks = []
    curr_chunk = []
    j = None
    d = None
    for (i,c) in zip(resnum,chain):
        if j is not None and \
           d is not None and \
           ( (i-1) != j or d != c ):
            chunks.append(curr_chunk)
            curr_chunk = []
        j = i
        d = c
        curr_chunk.append( (i,c) )
    chunks.append(curr_chunk)
    n_jumps = len(chunks) - 1
    if n_jumps > 0:
        for i in xrange(n_jumps):
            #obligate pairs
            new_pos1 = min( chunks[i][-1], chunks[i+1][0] )
            new_pos2 = max( chunks[i][-1], chunks[i+1][0] )
            new_pair = [new_pos1,new_pos2]
            if not already_listed_in_obligate_pair( new_pair, obligate_pair ):
                obligate_pair.extend( new_pair )

#######################################################################
working_input_res = []
working_input_chain = []
for m,chain in zip(input_res,input_chain):
    if (m,chain) not in zip(working_res,working_chain):
        # ignore chain
        if all( map( lambda x: x == '', working_chain ) ): chain = ''
    if (m,chain) not in zip(working_res,working_chain):
        raise ValueError('Input residue %s,%s not in working_res!!' % (m,chain) )
    i = zip(working_res,working_chain).index( (m,chain) )
    working_input_res.append( i+1 )
    working_input_chain.append( chain )

canonical_pairs = flatten( get_all_stems(working_secstruct) )
general_pairs = flatten( get_all_stems(working_secstruct_general))

for p in general_pairs:
    if p not in canonical_pairs:
        new_pair = [ (working_res[p[0]-1],working_chain[p[0]-1]),
                     (working_res[p[1]-1],working_chain[p[1]-1]) ]
        if not already_listed_in_obligate_pair( new_pair, obligate_pair ): obligate_pair.extend( new_pair )

fasta_file_outstring, params_file_outstring = \
    prepare_fasta_and_params_file_from_sequence_and_secstruct(
    working_sequence, working_secstruct, fixed_stems, working_input_res )

working_cst_file = ""
cst_file_outstring = ""
if len( cst_file ) > 0: working_cst_file = tag+"_"+cst_file

# also create a constraints file that will help close chains across 2-5 residue gaps...
# this is a slight hack, but if it works, might be worth putting a term into Rosetta, as
# well as automated handling of "working_res"
cst_gaps = []
for m in range( len(working_res)-1 ):
    gap = working_res[m+1] - working_res[m]
    if ( gap > 1 and gap < 6 ): cst_gaps.append( [ m+1, m+2, gap, working_res[m], working_res[m+1] ] )
#print "GAPS TO APPLY CST: ", cst_gaps

if len( cst_gaps ) > 0 and is_cst_gap:

    if len( working_cst_file ) == 0:  working_cst_file = tag+'_closegaps.cst'
    cst_file_outstring += "[ atompairs ]\n"
    for cst_gap in cst_gaps:
        stdev = 10.0
        gap_length = cst_gap[2]
        max_dist = gap_length * 5.0 + 4
        bonus = 200.0
        cst_file_outstring +=  " O3* %d  C5* %d   FADE %6.3f  %6.3f  %6.3f %6.3f %6.3f \n" % \
            ( cst_gap[0], cst_gap[1],  -stdev, max_dist, stdev, -1*bonus, bonus)


working_data_file = ""
if len( data_file ) > 0:
    lines = open( data_file ).readlines()
    working_data_string = ""
    for line in lines:
        cols = string.split( line.replace( '\n','') )
        if cols[0] == 'DMS':
            data_res = int( cols[1] )
            working_data_res = working_res_map( [data_res], working_res )
            if len( working_data_res ) > 0: working_data_string += cols[0]+' '+make_tag(working_data_res)+' '+string.join(cols[2:]) + '\n'
        else:
            data_res = map( lambda x:int(x), cols[1:] )
            working_data_res = working_res_map( data_res, working_res )
            if len( working_data_res ) > 0: working_data_string += cols[0]+' '+make_tag(working_data_res)+'\n'
    if len( working_data_string ) > 0:
        working_data_file = tag + "_" + data_file
        fid = open( working_data_file, 'w' )
        fid.write( working_data_string )
        fid.close()


assert( is_even( len(obligate_pair) ) )
if obligate_pair and len( obligate_pair[ 0 ] ) > 0:
    working_res_chain = zip( working_res, working_chain )
    for m in range( len( obligate_pair)/2 ):
        pos1 = obligate_pair[ 2*m ]
        pos2 = obligate_pair[ 2*m+1 ]
        if pos1 not in working_res_chain: continue
        if pos2 not in working_res_chain: continue
        working_pos1 = working_res_chain.index( pos1 ) + 1
        working_pos2 = working_res_chain.index( pos2 ) + 1
        params_file_outstring += "OBLIGATE  PAIR %d %d X X X \n" % (working_pos1, working_pos2)


assert( is_even( len(remove_obligate_pair) ) )
if len( remove_obligate_pair ) > 0:
    for m in range( len( remove_obligate_pair)/2 ):
        pos1 = remove_obligate_pair[ 2*m ]
        pos2 = remove_obligate_pair[ 2*m+1 ]
        if pos1 not in working_res: continue
        if pos2 not in working_res: continue
        pos1 = working_res.index( pos1 ) + 1
        pos2 = working_res.index( pos2 ) + 1
        params_file_lines = params_file_outstring.split( '\n' )
        new_lines = []
        for line in params_file_lines:
            found_pair = 0
            if line.find( 'OBLIGATE' ) > -1:
                pairs = line.split( 'PAIR' )[1:]
                for pair in pairs:
                    cols = pair.split()
                    if int(cols[0]) == pos1 and int(cols[1]) == pos2:
                        found_pair = 1
                        break
                    if int(cols[0]) == pos2 and int(cols[1]) == pos1:
                        found_pair = 1
                        break
            if found_pair:
                print 'Removing line: ', line
                continue
            new_lines.append( line )
        params_file_outstring = string.join( new_lines, '\n' )


assert( ( len(obligate_pair_explicit) % 5 == 0 ) )
if len( obligate_pair_explicit ) > 0:
    for m in range( len( obligate_pair_explicit)/5 ):
        pos1 = int( obligate_pair_explicit[ 5*m ] )
        pos2 = int( obligate_pair_explicit[ 5*m+1 ] )
        edge1 = obligate_pair_explicit[ 5*m+2 ]
        edge2 = obligate_pair_explicit[ 5*m+3 ]
        orientation = obligate_pair_explicit[ 5*m+4 ]
        if pos1 not in working_res: continue
        if pos2 not in working_res: continue
        pos1 = working_res.index( pos1 ) + 1
        pos2 = working_res.index( pos2 ) + 1
        params_file_outstring += "OBLIGATE  PAIR %d %d %s %s %s \n" % (pos1, pos2, edge1, edge2, orientation)


if len( chain_connection ) > 0:
    assert( len( chain_connection ) >= 4 )
    if 'SET1' in chain_connection: # new format is more flexible.
        which_set = 0
        chain_connection_sets = []
        resnum1 = []
        resnum2 = []
        for k in range( len( chain_connection) + 1 ):
            if k == len( chain_connection ) or chain_connection[k] == 'SET1':
                working_resnum1 = working_res_map( resnum1, working_res )
                working_resnum2 = working_res_map( resnum2, working_res )
                if len( working_resnum1 ) > 0 and len( working_resnum2 ) > 0: chain_connection_sets.append( [working_resnum1, working_resnum2] )
                resnum1 = []
                resnum2 = []
                which_set = 1
                continue
            if chain_connection[k] == 'SET2':
                which_set = 2
                continue
            assert( which_set > 0 )
            if which_set == 1: get_ints( chain_connection[k], resnum1 )
            if which_set == 2: get_ints( chain_connection[k], resnum2 )
        if len( chain_connection_sets ) > 0:
            params_file_outstring += "CHAIN_CONNECTION"
            for chain_connection_set in chain_connection_sets: params_file_outstring += "   SET1%s SET2%s" % (make_tag_with_dashes(chain_connection_set[0]), make_tag_with_dashes(chain_connection_set[1]) )
            params_file_outstring += "\n"
    else: # legacy format
        chain_connection = map( lambda x:int(x), chain_connection )
        assert( len( chain_connection ) % 4 == 0 )
        n_connect = len( chain_connection ) / 4
        working_chain_connection = working_res_map( chain_connection, working_res )
        for i in xrange(n_connect):
            curr_0 = i * 4
            seg1_start = working_chain_connection[curr_0]
            seg1_stop  = working_chain_connection[curr_0 + 1]
            seg2_start = working_chain_connection[curr_0 + 2]
            seg2_stop = working_chain_connection[curr_0 + 3]
            params_file_outstring += "CHAIN_CONNECTION SEGMENT1 %d %d  SEGMENT2 %d %d \n" % (seg1_start, seg1_stop, seg2_start, seg2_stop )

# need to handle Mg(2+)
#mg_pos = []
#for i in range( len( sequence ) ):
#    if sequence[i]=='z': mg_pos.append( i+1 )
#working_mg_pos = working_res_map( mg_pos, working_res )
#if len( working_mg_pos ) > 0:
#    for i in mg_pos:
#        cutpoint_open_res_chain.append( i-1 )
#        virtual_anchor.append( i )
def get_rid_of_previously_defined_cutpoints( cutpoint, working_res, working_chain ):
    cutpoint_filter = []
    for m in cutpoint:
        if ( m == len( working_res ) ): continue
        if ( working_res[ m-1 ]+1 != working_res[ m ] ): continue # already recognized as a cutpoint.
        if ( len( working_chain ) > 0 and working_chain[ m-1 ] != working_chain[ m ] ): continue # already recognized as a cutpoint.
        cutpoint_filter.append( m )
    return cutpoint_filter

if len( cutpoint_open_res_chain[0] ) > 0:
    working_cutpoint_open = working_res_map( cutpoint_open_res_chain, working_res, working_chain )
    working_cutpoint_open = get_rid_of_previously_defined_cutpoints( working_cutpoint_open, working_res, working_chain )
    if len( working_cutpoint_open ) > 0:    params_file_outstring += "CUTPOINT_OPEN  "+make_tag( working_cutpoint_open )+ "\n"

if len( cutpoint_closed_res_chain[0] ) > 0:
    working_cutpoint_closed = working_res_map( cutpoint_closed_res_chain, working_res, working_chain )
    working_cutpoint_closed = get_rid_of_previously_defined_cutpoints( working_cutpoint_closed, working_res, working_chain )
    if len( working_cutpoint_closed ) > 0:    params_file_outstring += "CUTPOINT_CLOSED  "+make_tag( working_cutpoint_closed )+ "\n"

if len( virtual_anchor_res_chain[0] ) > 0:
    working_virtual_anchor = working_res_map( virtual_anchor_res_chain, working_res, working_chain )
    params_file_outstring += "VIRTUAL_ANCHOR "+make_tag( working_virtual_anchor )+ "\n"


if len(cst_file) > 0:  # also have input data...
    lines = open( cst_file ).readlines()

    for line in lines:
        if len( line ) == 0: continue
        if line[0] == "[":
            cst_tag = line.split()[1]
            if ( cst_tag == "coordinates" ):
                cst_file_outstring += line
            else:
                assert( cst_tag == "atompairs" ) # must have cst file with atompairs first, then coordinates.
                if len( cst_file_outstring ) == 0:
                    cst_file_outstring += line
                print "WARNING! WARNING! WARNING! "
                print "if you have atompairs in your .cst file, make sure to put them first!"
        else:
            cols = line.split()
            if len( cols ) < 4: continue
            if int(cols[1]) not in working_res: continue
            if int(cols[3]) not in working_res: continue
            newline = "%s %d %s %d  %s\n" % (cols[0],  working_res.index(int(cols[1]))+1, cols[2],  working_res.index( int(cols[3]) )+1, string.join( cols[4:] ) )
            cst_file_outstring += newline

fasta_file = tag+'.fasta'
print
print 'Writing to fasta file: ', fasta_file
fid = open( fasta_file, 'w' )
#print fasta_file_outstring
fid.write( fasta_file_outstring )
fid.close()

params_file = tag+'.params'
print 'Writing to params file: ', params_file
fid = open( params_file, 'w' )
#print params_file_outstring
fid.write( params_file_outstring )
fid.close()

if len( working_cst_file ) > 0 :
    cst_file = tag+'.cst'
    print 'Writing to cst file: ', working_cst_file
    fid = open( working_cst_file, 'w' )
    #print cst_file_outstring
    fid.write( cst_file_outstring )
    fid.close()

assert( not ( len( native_pdb )>0 and len( working_native_pdb ) > 0 ) )
if ( len(native_pdb) > 0 and len( working_res ) > 0):
    assert( exists( native_pdb ) )
    command = "pdbslice.py " + native_pdb + " -subset"
    for m in working_res: command += " %d" % m
    command += " "+tag+"_"
    system( command )
    working_native_pdb = "%s_%s" % (tag,native_pdb)
    print "Writing native to:", working_native_pdb


#########################################
print
print "Sample command line: "


if rosetta_folder == "":
	rosetta_folder = None
if extension == "":
	extension = None

command  = rosetta_exe('rna_denovo', rosetta_folder, extension)
command += " -nstruct %d -params_file %s -fasta %s  -out:file:silent %s.out -include_neighbor_base_stacks " % (nstruct, params_file, fasta_file, tag )
if no_minimize:
    command += " -minimize_rna false"
else:
    command += " -minimize_rna true"

if len( working_native_pdb ) > 0:
    command += " -native %s " % working_native_pdb
elif len( native_pdb ) > 0:
    command += " -native %s " % native_pdb

if len( extra_minimize_res_chain[0] ) > 0:
    working_extra_minimize_res = working_res_map( extra_minimize_res_chain, working_res, working_chain )
    command += " -extra_minimize_res " + make_tag_with_dashes( working_extra_minimize_res, [] )

if len( extra_minimize_chi_res_chain[0] ) > 0:
    working_extra_minimize_chi_res = working_res_map( extra_minimize_chi_res_chain, working_res, working_chain )
    command += " -extra_minimize_chi_res " + make_tag_with_dashes( working_extra_minimize_chi_res, [] )

if len( input_pdbs ) > 0:
    command += " -s"
    for pdb in input_pdbs:
        assert( exists( pdb ) )
        command += " "+pdb
    assert( len( input_res ) > 0 )

if len( input_silent_files ) > 0:
    command += " -in:file:silent"
    for silent_file in input_silent_files:
        assert( exists( silent_file ) )
        command += " "+silent_file
    assert( len( input_res ) > 0 )

if len( input_res ) > 0:
    command += " -input_res " +make_tag_with_dashes( working_input_res, [] )

if len( working_cst_file ) > 0:
    command += " -cst_file " + working_cst_file

if len( working_data_file ) > 0:
    command += " -data_file " + working_data_file

if bps_moves:
    command += " -bps_moves "

command += ' ' + extra_args

command += ' -output_res_num ' + make_tag_with_dashes( working_res, working_chain )

print command

print "outputting command line to: ", out_script
with open( out_script, 'w' ) as fid:
    fid.write( command + "\n" )

print
print
print "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!"
print " This script rna_denovo_setup.py is deprecated.  "
print " Instead, run your flags straight into the rosetta rna_denovo application: "
print
print "    rna_denovo",
for arg in argv_input[1:]:
    if arg == '-no_minimize': continue
    print arg,
command = ""
command += " -include_neighbor_base_stacks "
if no_minimize:
    command += " -minimize_rna false"
else:
    command += " -minimize_rna true"
print command
print
print "WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING! WARNING!"


