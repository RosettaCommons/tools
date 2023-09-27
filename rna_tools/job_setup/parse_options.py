#!/usr/bin/env python
import string

def get_resnum_chain( input_string, resnums, chains, segids ): # could be of the form A:1-4 or A1-4 or 1-4

    if ( input_string.find( ',' ) > -1 ):
        for substring in input_string.split(','): get_resnum_chain( substring, resnums, chains, segids )
        return True
    #assert( len( input_string ) > 0 )
    #if len(segids) == 0: segids = ['    ' for _ in resnums]
    assert( len( resnums ) == len( chains ) )
    assert( len( resnums ) == len( segids ) )
    ok = True

    int_string = input_string
    chain = ''
    segid = '    '

    # extract chain info, if present.
    # colon is clear indicator of chain
    if input_string.count(':') == 2:
        first = input_string.find( ':' )
        chain = input_string[:first]
        second = first+2+input_string[ first+2 : ].find(':')
        segid = input_string[ first+1 : second ]
        int_string = input_string[second+1: ]
    else:
        colon_pos = input_string.find( ':' )
        if colon_pos > -1:
            chain = input_string[:colon_pos]
            int_string = input_string[ colon_pos+1 : ]
        else: # look for non-numeric characters
            pos = 0
            while pos < len( input_string ):
                try:
                    if ( input_string[pos] != '-' ): blah = int( input_string[ pos ] )
                    break
                except:
                    chain += input_string[ pos ]
                pos = pos + 1
            int_string = input_string[ pos: ]

    ints = []
    if len( int_string ) == 0 and len( chain ) == 1: # special case, get the whole chain.
        resnums.append( 'all' ) # means get everything from this chain
        chains.append( chain )
    else:
        ok = get_ints( int_string, ints )
        for i in range( len( ints ) ):
            resnums.append( ints[i] )
            chains.append( chain )
            segids.append( segid )
    assert( len( resnums ) == len( chains ) )
    assert( len( resnums ) == len( segids ) )
    return ok


def get_ints( int_string, value ): # could be of the form 5-8  or -9--5
    if ( len( int_string ) == 0 ): return False

    # find hyphen separator
    dashpos = 0
    if int_string[0] == '-': dashpos = 1  # can't be the first dash
    while dashpos < len( int_string ) and int_string[ dashpos ] != '-': dashpos += 1
    if dashpos == len( int_string ):
        try:
            value.append( int( int_string ) )
        except:
            return False
    else:
        try:
            int_start = int( int_string[          :dashpos] )
            int_stop  = int( int_string[dashpos+1 :       ] )
            assert( int_stop > int_start )
            for m in range( int_start, int_stop+1): value.append( m )
        except:
            return False
    return True

def has_repeated_flags( argv ):
    flags = []
    for arg in argv:
        if arg[0] == "-":
            if arg in flags: return True
            flags.append( arg )
    return False

################################################################
def parse_options( argv, tag, default):

    assert( not has_repeated_flags( argv ) )

    value = default

    if isinstance( default, list ): value = []
    RESNUM_CHAIN_SEGID = isinstance( default, list ) and len( default ) == 3 and \
                   isinstance( default[0], list ) and \
                   isinstance( default[0][0], int ) and \
                   isinstance( default[1], list ) and \
                   isinstance( default[1][0], str )
    if RESNUM_CHAIN_SEGID: value = [ [], [], [] ]; # need three lists.

    if argv.count( "-"+tag ):

        pos = argv.index( "-"+tag )
        del( argv[ pos ] )

        if ( ( default == 0 or isinstance( default, bool ) ) and
             ( pos == (len( argv ) ) or
               argv[ pos ][0] == '-' ) ): # Just a boolean
            value = 1
        elif( isinstance( default, list ) ):
            while ( pos < len (argv) and not (argv[ pos ][0] == '-') ):
                if isinstance( default[0], int ):
                    ok_int_string = get_ints( argv[pos], value )
                    if not ok_int_string: break
                elif isinstance( default[0], float ):
                    value.append( float( argv[ pos ] ) )
                elif RESNUM_CHAIN_SEGID:
                    ok = get_resnum_chain( argv[pos], value[0],value[1],value[2] )
                    if not ok: break
                else:
                    value.append( argv[ pos ] )
                del( argv[ pos ] )
        elif isinstance( default, int ):
            value = int( argv[ pos ] )
            del( argv[ pos ] )
        elif isinstance( default, float ):
            value = float( argv[ pos ] )
            del( argv[ pos ] )
        else:
            value = argv[ pos ]
            del( argv[ pos ] )

    #print tag, value

    return value

