#!/usr/bin/python

################################################################
def parse_options( argv, tag, default):
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
                    value.append( int( argv[ pos ] ) )
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

