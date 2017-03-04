#!/usr/bin/env python

from SWA_util import *
################################################################
def parse_options( argv, tag, default):

	print 'FOUND: ', tag


	if argv.count( "-"+tag ):  ###Found the option

		actual_offset=1

		pos = argv.index( "-"+tag )   ###Position of the option name

		if ( ( pos == (len( argv )-1) or argv[ pos+1 ][0] == '-' ) ):
			error_exit_with_message("Invalid parse_option input")
		elif (default=="False" or default=="True"): # Python boolean
			print "%s=%s" %(tag, argv[ pos + 1 ])
			if ( argv[ pos + 1 ] == "True"):
				value=True
			elif ( argv[ pos + 1 ] == "False"):
				value=False
			else:
				error_exit_with_message('(%s != "True") and  (%s != "False")' % (tag, tag))
		elif (default=="false" or default=="true"): # C++ boolean string
			print "%s=%s" %(tag, argv[ pos + 1 ])
			if ( argv[ pos + 1 ] == "true"):
				value="true"
			elif ( argv[ pos + 1 ] == "false"):
				value="false"
			else:
				error_exit_with_message('(%s != "true") and  (%s != "false")' % (tag, tag))
		elif( isinstance( default, str ) ): #normal string
			print "%s=%s" %(tag, argv[ pos + 1 ])
			value=argv[ pos + 1 ]

		elif( isinstance( default, list ) ):
			value = []
			offset = 1
			while ( (pos + offset) < len (argv) and not (argv[ pos+offset ][0] == '-') ):
				actual_offset=offset
				if isinstance( default[0], int ):
					value.append( int( argv[ pos+offset ] ) )
				elif isinstance( default[0], float ):
					value.append( float( argv[ pos+offset ] ) )
				else:
					value.append( argv[ pos+offset ] )
				offset += 1
		elif isinstance( default, int ):
			value = int( argv[ pos + 1 ] )
		elif isinstance( default, float ):
			value = float( argv[ pos + 1 ] )
		else:
			value = argv[ pos + 1 ]

		for index in range(pos+actual_offset, pos-1, -1):  #index from pos+actual_offset to pos
			del(argv[index])

		return value

	else: #Return the default value
		if(default=="False" or default=="True"): # Python boolean
			print "%s=%s" %(tag, default)
			if( default == "True"):
				return True
			elif( default == "False"):
				return False
			else:
				error_exit_with_message('(%s != "True") and  (%s != "False")' % (tag, tag))
		elif( isinstance( default, list ) and isinstance( default[0], int ) ):
			default = [] #A empty list, why do we need this??
			return default
		else:
			return default




