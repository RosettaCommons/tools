# -*- coding: utf-8 -*-
# :noTabs=true:

import string

from os import popen,system
from os.path import exists,dirname,basename,abspath

#from rna_conversion import make_rna_rosetta_ready
from get_sequence import get_sequence


MAX_SEQUENCE_LENGTH = 32

nts = ['g','c','u','a','G','C','U','A','z','Z']
secstruct_chars = ['(',')','[',']','{','}','.']
spacers = ['+','*',' ',','] # any of these are OK as strand separators
complement = {'a':['u'], 'u':['a','g'], 'c':['g'], 'g':['c','u']};

def ValidationError( string ):
    print string
    exit()

def join_sequence( sequence ):
    sequence_joined = ''
    chainbreak_pos = []
    count = 0
    for m in range( len( sequence ) ):
        c = sequence[m]
        if c in nts or c in secstruct_chars:
            sequence_joined += c.lower()
            count += 1
        elif c in spacers:
            if ( c == 0 ):
                raise ValidationError( "Cannot start secstruct with spacer!" )
                return None
            chainbreak_pos.append( count )
        else:
            raise ValidationError( "Unrecognized character in sequence: %s" % c  )
            return None
    return ( sequence_joined, chainbreak_pos )

def prepare_fasta_and_params_file_from_sequence_and_secstruct( sequence, secstruct='', fixed_stems = False, input_res = None ):

    fasta_file_outstring  = ""
    params_file_outstring = ""

    if len( secstruct ) > 0 and len( secstruct ) != len( sequence ):
        raise ValidationError("If secstruct is specified, length of secstruct must be equal to length of sequence.!" )
        return None

    if len(sequence.split('\n')) != 1:
        raise ValidationError("sequence must be one line!")
        return None

    if len( secstruct ) > 0 and len(secstruct.split('\n')) != 1:
       #raise ValidationError("secstruct must be blank or one line!")
        return None

    # find chainbreaks
    ( sequence_for_fasta, chainbreak_pos ) = join_sequence( sequence )

    fasta_file_outstring = "> Input sequence: " + sequence + "\n"
    fasta_file_outstring += sequence_for_fasta
    fasta_file_outstring = convert_fasta_to_rosetta_format( fasta_file_outstring )

    if len( chainbreak_pos ) > 0  and len( secstruct ) == 0:
        raise ValidationError( "If multiple strands, must specify some kind of Watson/Crick pairing to connect the strands in secstruct. [Alternatively, advanced users can set up a params file to input non-Watson/Crick pairing to connect strands.]!" )
        return None

    # Figure out stems & cutpoints for params file
    if len( secstruct ) > 0:
        params_file_outstring += "# Input secstruct " + secstruct + "\n"

        if len( chainbreak_pos ) > 0:
            params_file_outstring += "CUTPOINT_OPEN "
            for m in chainbreak_pos: params_file_outstring += " %d" % m
            params_file_outstring += "\n"

        #Changes to meaning of 'fixed_stems' here
        stems = get_stems( secstruct, chainbreak_pos, '(', ')', sequence_for_fasta )
        params_file_outstring += output_stems( stems, fixed_stems, input_res )

        stems = get_stems( secstruct, chainbreak_pos, '[', ']', sequence_for_fasta )
        params_file_outstring += output_stems( stems, fixed_stems, input_res )

        stems = get_stems( secstruct, chainbreak_pos, '{', '}', sequence_for_fasta )
        params_file_outstring += output_stems( stems, fixed_stems, input_res )

    #assume that 'z' means Mg(2+)
    fasta_file_outstring = fasta_file_outstring.replace( 'z','Z[MG]')

    return ( fasta_file_outstring, params_file_outstring )


def output_stems( stems, fixed_stems = False, input_res = None ):
    outstring = ''
    for stem_res in stems:
        outstring += 'STEM '
        for pair in stem_res:
             outstring += ' PAIR %d %d W W A ' % (pair[0], pair[1])
        outstring += '\n'
        if fixed_stems and input_res is not None:
            for pair in stem_res:
                if pair[0] in input_res and pair[1] in input_res:
                    break
            else:
                pair = stem_res[0]
                outstring += 'OBLIGATE PAIR %d %d W W A \n' % (pair[0], pair[1])
    return outstring

def get_all_stems( secstruct, chainbreak_pos = [], sequence_for_fasta='' ):
    #chainbreak_pos = []
    #sequence_for_fasta = ''
    stems = []
    for stem in get_stems( secstruct, chainbreak_pos, '(', ')', sequence_for_fasta ): stems.append( stem )
    for stem in get_stems( secstruct, chainbreak_pos, '[', ']', sequence_for_fasta ): stems.append( stem )
    for stem in get_stems( secstruct, chainbreak_pos, '{', '}', sequence_for_fasta ): stems.append( stem )
    return stems

def get_stems( line, chainbreak_pos, left_bracket_char = '(', right_bracket_char = ')', sequence_for_fasta='' ):
    count = 0
    left_brackets = []
    pair_map = {}
    all_pairs = []

    for i in range( len(line) ):
        if line[i] in spacers: continue
        count += 1
        if line[i] == left_bracket_char:  left_brackets.append( count )
        if line[i] == right_bracket_char:
            if len( left_brackets ) == 0:
                raise ValidationError( "Number of right brackets does not match left brackets" )
            res1 = left_brackets[-1]
            res2 = count
            del( left_brackets[-1] )
            pair_map[ res1 ] = res2
            pair_map[ res2 ] = res1
            all_pairs.append( [res1,res2] )
            if len( sequence_for_fasta ) > 0 and not ( sequence_for_fasta[res1-1] in complement[ sequence_for_fasta[res2-1] ] ):
                #raise ValidationError( "Not complementary at positions %s%d and %s%d!"  % (sequence_for_fasta[res1-1],res1,sequence_for_fasta[res2-1],res2) )
                print "[WARNING] Not complementary at positions %s%d and %s%d!"  % (sequence_for_fasta[res1-1],res1,sequence_for_fasta[res2-1],res2) )
                continue

    if ( len (left_brackets) > 0 ):
        raise ValidationError( "Number of right brackets does not match left brackets" )
    numres = count

    # Parse out stems
    already_in_stem = {}
    for i in range( numres ): already_in_stem[ i ] = 0

    stems = []
    stem_count = 0
    for i in range( 1, numres+1 ):
        if pair_map.has_key( i ) and not already_in_stem[ i ]:  # In a base pair
            k = i
            stem_count += 1
            stem_res = []

            stem_res.append( [k, pair_map[k]] )
            already_in_stem[ k ] = 1
            already_in_stem[ pair_map[k] ] = 1

            # Can we extend in one direction?
            while( pair_map.has_key( k + 1 ) and  pair_map[ k+1 ] == pair_map[ k ] - 1  and not already_in_stem[k+1] and (not k in chainbreak_pos) and ( not pair_map[k+1] in chainbreak_pos )  ):
                k += 1
                stem_res.append( [k, pair_map[k]] )
                already_in_stem[ k ] = 1
                already_in_stem[ pair_map[k] ] = 1

            stems.append( stem_res )
    return stems

def convert_fasta_to_rosetta_format( fasta_file_sequence ):
    '''
    # read lines from fasta file
    if not exists( fasta_file_name ):
        stderr.write( 'Cannot find filename: %s.\n' % fasta_file_name )
        return None

    lines = open( fasta_file_name ).readlines()
        '''

    lines = map(lambda x: x+'\n', fasta_file_sequence.split('\n') )

    line = lines[0]
    if not line[0] == '>':
        #stderr.write( 'First line of fasta must start with \'>\'\n' )
        raise ValidationError("First line of fasta must start with '>'!")
        return None

    sequence = string.join( lines[1:] ).replace( ' ','').replace('\n','')
    sequence = sequence.lower()

    if len( sequence ) > MAX_SEQUENCE_LENGTH:
        #stderr.write( 'Cannot model a sequence with length greater than %d. Your sequence has length %d.\n' % (MAX_SEQUENCE_LENGTH, len( sequence) ) )
        #raise ValidationError( 'Cannot model a sequence with length greater than %d. Your sequence has length %d.!' % (MAX_SEQUENCE_LENGTH, len( sequence) ) )
        #return None
        print ''

    # quick validation
    goodchars = ['a','c','g','u','z']
    for char in sequence:
        if char not in goodchars:
            #stderr.write( 'Character %s is not a,c,g, or u\n'% char )
            raise ValidationError( 'Character %s is not a,c,g, or u!'% char )
            return None

    # here's the output
    outstring = ''
    outstring += line
    outstring += sequence
    outstring += '\n'

    return outstring


def convert_pdb_to_rosetta_format( pdb_file_name ):
    outstring = make_rna_rosetta_ready( pdb_file_name )
    return outstring


def does_PDB_match_fasta(pdb, fasta_file_sequence ):
    '''
    if not exists( PDB_file_name  ):
        stderr.write( 'Could not find file name: %s.\n' % PDB_file_name )
        return None

    if not exists( fasta_file_name  ):
        stderr.write( 'Could not find file name: %s.\n' % fasta_file_name )
        return None

    seq_PDB = get_sequence( PDB_file_name )
    seq_fasta = string.join( open( fasta_file_name ).readlines()[1:] ).replace( '\n', '' ).replace(' ','' )
    return ( seq_PDB == seq_fasta ) '''

    seq_PDB = get_sequence(pdb)
    seq_fasta = string.join( fasta_file_sequence.split('\n')[1:] ).replace( '\n', '' ).replace(' ','' )

    #print 'does_PDB_match_fasta: pdb=%s\nfasta=%s' % (pdb, fasta_file_sequence)
    #print 'does_PDB_match_fasta: %s %s' % (seq_PDB, seq_fasta)

    return ( seq_PDB == seq_fasta )


longer_names={'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
              'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
              'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
              'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
              'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
              ' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u' }


# accepts a pdb file name, returns a string with pdb entries -- or None if there is an error.
def make_rna_rosetta_ready( pdb, removechain=False, ignore_chain=True, chainids = [], no_renumber = False, removeions = False, old_rna = False ):

    #fastaid = stderr
    num_model = 0
    #max_model = 60 # for virus
    max_model = 0 # for virus

    outstring = ''

    lines = map(lambda x: x+'\n', pdb.split('\n') )

    #fastaid.write('>'+pdbname+'\n');

    oldresnum = '   '
    count = 0;

    for i in range( len( chainids ) ) :
        if chainids[i] == '_':
            chainids[i] = ' '

    goodnames = ['  A','  C','  G','  U',' rA',' rC',' rG',' rU',' MG', ' IC',' IG']
    hetatm_map = { '5BU':'  U', ' MG':' MG', 'OMC':'  C', '5MC':'  C', 'CCC':'  C', ' DC':'  C', 'CBR':'  C', 'CBV':'  C', 'CB2':'  C', '2MG':'  G', 'H2U':'  U', 'PSU':'  U', '5MU':'  U', 'OMG':'  G', '7MG':'  G', '1MG':'  G', 'GTP':'  G', 'AMP':'  A', ' YG':'  G', '1MA':'  A', 'M2G':'  G', 'YYG':'  G', ' DG':'  G', 'G46':'  G', ' IC':' IC',' IG':' IG' }


    if removeions:  goodnames.remove(' MG')


    for line in lines:
        if len(line)>5 and line[:6]=='ENDMDL':
            num_model += 1
            if num_model > max_model:  break #Its an NMR model.
        if len(line) <= 21:  continue
        if (line[21] in chainids or ignore_chain):
            line_edit = line

            if line[0:3] == 'TER' and False:
                continue
            elif (line[0:6] == 'HETATM') & (line[17:20]=='MSE'): #Selenomethionine
                line_edit = 'ATOM  '+line[6:17]+'MET'+line[20:]
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit)>75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]
            elif (line[0:6] == 'HETATM') & ( line[17:20] in hetatm_map.keys()):
                line_edit = 'ATOM  '+line[6:17]+ hetatm_map[line[17:20]] + line[20:]

            #Don't save alternative conformations.
            if line[16] == 'A':
                continue;

            if line_edit[0:4] == 'ATOM':
                if line_edit[22] == ' ': resnum = line_edit[23:26]
                else:                    resnum = line_edit[22:26]
                if not resnum == oldresnum: #  or line_edit[12:16] == ' P  ':
                    longname = line_edit[17:20]
                    if longname == 'GTP':
                        longname = '  G'
                    elif longname == 'G  ':
                        longname =   '  G'
                    elif longname == ' DG':
                        longname =   '  G'
                    elif longname == 'A  ':
                        longname =   '  A'
                    elif longname == 'C  ':
                        longname =   '  C'
                    elif longname == 'U  ':
                        longname =   '  U'
                    elif longname == 'GUA':
                        longname = '  G'
                    elif longname == 'ADE':
                        longname = '  A'
                    elif longname == 'CYT':
                        longname = '  C'
                    elif longname == 'URA':
                        longname = '  U'
                    elif longname == 'URI':
                        longname = '  U'
                    elif longname == ' rA':
                        longname =   '  A'
                    elif longname == ' rC':
                        longname =   '  C'
                    elif longname == ' rU':
                        longname =   '  U'
                    elif longname == ' rG':
                        longname =   '  G'
                    elif longname == '  I':
                        longname =   '  G'
                    else:
                        if longname not in goodnames:    continue

                    if longer_names.has_key(longname):
                        #fastaid.write( longer_names[longname] );
                        pass
                    else:
                        #fastaid.write( 'X')
                        pass

                    #print "AAH ==> " ,  resnum, oldresnum, line_edit
                    count = count + 1

                oldresnum = resnum

                if not longname in goodnames:
                    continue

                newnum = '%4d' % count
                if no_renumber: newnum = '%4s' % resnum

                line_edit = line_edit[0:16] + ' ' + longname + line_edit[20:22] + \
                            newnum + line_edit[26:]
                if removechain:
                    line_edit = line_edit[0:21]+'  '+line_edit[23:]

                line_edit = line_edit.replace('2HO*', "HO2'")
                line_edit = line_edit.replace('5HO*', "HO5'")
                line_edit = line_edit.replace('2H5*', "H5''")

                line_edit = line_edit.replace('*', "'")
                line_edit = line_edit.replace('O1P', 'OP1')
                line_edit = line_edit.replace('O2P', 'OP2')


                if old_rna:
                    line_edit = line_edit.replace('  A', ' rA')
                    line_edit = line_edit.replace('  C', ' rC')
                    line_edit = line_edit.replace('  G', ' rG')
                    line_edit = line_edit.replace('  U', ' rU')

                outstring += line_edit

    #fastaid.write('\n')

    return outstring

