#!/usr/bin/python


def get_ints( int_string, value ): # could be of the form 5-8  or -9--5
    assert( len( int_string ) > 0 )

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

    if argv.count( "-"+tag ):

        pos = argv.index( "-"+tag )
        del( argv[ pos ] )

        if ( ( default == 0 or isinstance( default, bool ) ) and
             ( pos == (len( argv ) ) or
               argv[ pos ][0] == '-' ) ): # Just a boolean
            value = 1
        elif( isinstance( default, list ) ):
            value = []
            while ( pos < len (argv) and not (argv[ pos ][0] == '-') ):
                if isinstance( default[0], int ):
                    ok_int_string = get_ints( argv[pos], value )
                    if not ok_int_string: break
                elif isinstance( default[0], float ):
                    value.append( float( argv[ pos ] ) )
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
    else:
        if isinstance( default, list ): value = []

    #print tag, value

    return value

