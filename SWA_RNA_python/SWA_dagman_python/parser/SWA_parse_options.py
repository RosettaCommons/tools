#!/usr/bin/env python

import copy
################################################################


from SWA_dagman_python.utility.SWA_util import *

################################################################
def parse_options( argv, tag, default, Verbose=True):

	if(argv.count( "-"+tag )>1): error_exit_with_message("argv.count( %s )>1 " %("-"+tag) ) #June 13, 2011

	if(argv.count( "-"+tag )):  ###Found the option

		actual_offset=1

		pos = argv.index( "-"+tag )   ###Position of the option name

		if (default=="False" or default=="True" or isinstance( default, bool ) ): # Python boolean
			# flag is 'solo' --> activate to True
			if ( pos == (len( argv )-1) or argv[ pos+1 ][0] == '-' ):
				if(Verbose): print "%s=%s" % (tag, "True")
				value = True
				actual_offset -= 1
			else: # value is given explicitly as True or False.
				if(Verbose): print "%s=%s" %(tag, argv[ pos + 1 ])
				if ( argv[ pos + 1 ] == "True"):
					value=True
				elif ( argv[ pos + 1 ] == "False"):
					value=False
				else:
					error_exit_with_message('(%s != "True") and  (%s != "False")' % (tag, tag))
		elif ( ( pos == (len( argv )-1) or argv[ pos+1 ][0] == '-' ) ):
			error_exit_with_message("Invalid parse_option input")
		elif (default=="false" or default=="true"): # C++ boolean string
			if(Verbose): print "%s=%s" %(tag, argv[ pos + 1 ])
			if ( argv[ pos + 1 ] == "true"):
				value="true"
			elif ( argv[ pos + 1 ] == "false"):
				value="false"
			else:
				error_exit_with_message('(%s != "true") and  (%s != "false")' % (tag, tag))

		elif( isinstance( default, str ) ): #normal string
			if(Verbose): print "%s=%s" %(tag, argv[ pos + 1 ])
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
#		print "OUTSIDE %s is list, default= " %tag, default

		if(default=="False" or default=="True"): # Python boolean
			if(Verbose): print "%s=%s" %(tag, default)
			if( default == "True"):
				return True
			elif( default == "False"):
				return False
			else:
				error_exit_with_message('(%s != "True") and  (%s != "False")' % (tag, tag))
		elif( isinstance( default, list ) and isinstance( default[0], int ) ):
			default = [] #A empty list, why do we need this??
			if(Verbose): print "%s is list, default= " %tag, default

			return default
		else:
			return default


################################################################
##Find a specific arg_tag in a string of args and replace the value.
def replace_arg_value(args_string, arg_tag, new_value, allow_missing=False):

	args_list=args_string.split()

	if(args_list.count( "-"+arg_tag )>1): error_exit_with_message("args_list.count( %s )>1 " %("-"+arg_tag) )

	if(allow_missing==False and (args_list.count( "-"+arg_tag )<1) ): error_exit_with_message("allow_missing==False and args_list.count( %s )<1 " %("-"+arg_tag) )

	found_arg_tag=False
	insert_pos=-1
	if(args_list.count( "-"+arg_tag )>0):
		insert_pos = args_list.index( "-"+arg_tag )
		found_arg_tag=True

	old_value=parse_options(args_list, arg_tag, new_value, Verbose=False) #Ok a little hacky..remove the old value including the arg_tag.

	#check that args_tag is removed from the args_list!
	if(args_list.count( "-"+arg_tag )>0): error_exit_with_message("args_list.count( %s )>0 AFTER removal!" %("-"+arg_tag) )

	if(found_arg_tag==False):
		print "found_arg_tag (%s)==False" %(arg_tag)
		insert_pos=len(args_list)

	return_args_list=[]
	for n in range(0, insert_pos):
		return_args_list.append(args_list[n])

	return_args_list.append("-"+arg_tag)

	if(isinstance( new_value, list ) ):
		return_args_list.extend(new_value)
	else:
		return_args_list.append(new_value)

	for n in range(insert_pos, len(args_list)):
		if(found_arg_tag==False): error_exit_with_message("found_arg_tag==False") #if found_arg_tag==False, then should add to end of list!
		return_args_list.append(args_list[n])

	if(isinstance( new_value, list ) ):
		if(len(return_args_list)!=(len(args_list)+len(new_value)+1)): error_exit_with_message("len(return_args_list)!=(len(args_list)+len(new_value)+1)")
	else:
		if(len(return_args_list)!=(len(args_list)+2)): error_exit_with_message("len(return_args_list)!=(len(args_list)+2)")

	return_args_string=list_to_string(return_args_list)

	return return_args_string

################################################################

def parse_seq_num_list_option( argv, tag, sort_list=True ):

	return parse_segment_string_list( parse_options(argv, tag, [""] ), do_sort=sort_list )

##############################################################################################################################
def ensure_no_duplicate_options(option_list_string, extra_existing_option_list=[]): #To should work for both python and Rosetta (C++) commands

	if(isinstance( option_list_string, str )==False):
		print "\n\nERROR option_list_string=", option_list_string
		error_exit_with_message("option_list_string is not a str!")

	option_list=option_list_string.split()

	option_names_so_far=[]

	for n in range(len(option_list)):
		if(option_list[n][0]=="-"):

			option_name=option_list[n]

			if(option_names_so_far.count(option_name)>0):
				exit_message ="\n\nERROR option_list_string=%s \n" %(option_list_string)
				exit_message+="option_name (%s) already exist in the option_names_so_far ( %s )\n" %(option_name, list_to_string(option_names_so_far, " | ") )
				error_exit_with_message(exit_message)

			if(extra_existing_option_list.count(option_name[1:])>0):
				exit_message ="\n\nERROR option_list_string=%s \n" %(option_list_string)
				exit_message+="option_name (%s) is in extra_existing_option_list (%s)\n!" %(option_name, list_to_string(extra_existing_option_list) )
				error_exit_with_message(exit_message)

			option_names_so_far.append(option_name)

##############################################################################################################################
def option_name_exist(option_list_string, option_name):

	if(isinstance( option_list_string, str )==False):
		print "\n\nERROR option_list_string=", option_list_string
		error_exit_with_message("option_list_string is not a str!")

	option_list=option_list_string.split()

	for n in range(len(option_list)):
		if(option_list[n][0]!="-"): continue

		if(option_name==option_list[n][1:]): return True

	return False

##############################################################################################################################
def get_option_name_args_safe(option_list_string, option_name, defualt): #get the option_args without removing the args from the option_list!

	if(isinstance( option_list_string, str )==False):
		print "\n\nERROR option_list_string=", option_list_string
		error_exit_with_message("option_list_string is not a str!")

	if(option_name_exist(option_list_string, option_name)==False):
		error_exit_with_message("option_name (%s) doesn't exist in option_list_string (%s) " %(option_name, option_list_string) )

	option_list=option_list_string.split()

	args_value=parse_options(option_list, option_name, defualt, Verbose=False)

	return args_value




