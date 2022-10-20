#!/usr/bin/env python

def make_tag( int_vector ):
    tag = ''
    for m in int_vector: tag += ' %d' % m
    return tag


def make_tag_with_dashes( int_vector, char_vector = 0, segid_vector = 0 ):
    tag = ''
    if not isinstance( char_vector, list ) or len( char_vector ) == 0:
         char_vector = []
         for m in range( len( int_vector ) ): char_vector.append( '' )
    if not isinstance( segid_vector, list ) or len( segid_vector ) == 0:
        segid_vector = []
        for m in range( len( int_vector ) ): segid_vector.append( "    " )

    assert( len( char_vector ) == len( int_vector ) )
    assert( len( segid_vector ) == len( int_vector ) )

    start_res = int_vector[0]
    for i in range( 1, len(int_vector)+1 ):
        if i==len( int_vector)  or \
                int_vector[i] != int_vector[i-1]+1 or \
                char_vector[i] != char_vector[i-1] or \
                segid_vector[i] != segid_vector[i-1]:

            stop_res = int_vector[i-1]
            tag += ' '
            if len( char_vector[i-1] ) > 0:
                assert( len( char_vector[i-1] ) == 1 )
                tag += char_vector[i-1]+':'
            if segid_vector[i-1] != "    ":
                print ( segid_vector )
                assert( len( segid_vector[i-1] ) == 4 )
                tag += segid_vector[i-1].strip()+":"
            if stop_res > start_res:
                tag += '%d-%d' % (start_res, stop_res )
            else:
                tag += '%d' % (stop_res )

            if ( i < len( int_vector) ): start_res = int_vector[i]

    return tag


def make_tag_with_conventional_numbering( int_vector, resnums, chains, segids ):
    tag_resnums = []
    tag_chains  = []
    tag_segids  = []
    for m in int_vector:
        tag_resnums.append( resnums[ m - 1 ] )
        tag_chains.append(  chains [ m - 1 ] )
        tag_segids.append(  segids [ m - 1 ] )
    return make_tag_with_dashes( tag_resnums, tag_chains, tag_segids )


def make_tag_with_dashes_and_commas( resnums, chains = 0, segids = 0 ):
    tag  = make_tag_with_dashes( resnums, chains, segids )
    tag = tag.replace(' ',',')
    if tag[0] == ',' or tag[0] == ' ' or tag[0] == '': tag = tag[1:]
    return tag

def make_tag_from_list_of_int_ranges( list_of_int_ranges ):
    #######################################################
    ### >>> make_tag_from_list_of_int_ranges( ['A:4-5', 'B:6-10'] )
    ### '  4 5  6 7 8 9 10'
    #######################################################
    tag = ''
    for int_range in list_of_int_ranges:
        tag +=' '+make_tag_from_int_range( int_range )
    return tag 


def make_tag_from_int_range( int_range ):
    #####################################
    ### >>> make_tag_from_int_range( 'B:4-5' )
    ### ' 4 5'
    ### Could ignoring the Chain ID be dangerous???
    #####################################
    if ':' in int_range:
        chain = int_range.split(':')[0]
        int_range = int_range.split(':')[1]
    else:  
        int_range = int_range
    if '-' in int_range:
        first_idx = int(int_range.split('-')[0])
        last_idx = int(int_range.split('-')[1])
    else:
        first_idx = int(int_range)
        last_idx = int(int_range)
    int_vector = [ x for x in xrange(first_idx, last_idx+1) ]
    return make_tag( int_vector )
