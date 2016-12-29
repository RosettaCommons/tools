#!/usr/bin/env python

from SWA_dagman_python.parser.SWA_parse_options import parse_options, replace_arg_value
from SWA_dagman_python.utility.SWA_util import *
from SWA_dagman_python.utility.DAGMAN_util import *


#############################################################################################################

def get_sampler_post_process_generic_file(output_foldername, prepend_str, append_str):

		if(isinstance( prepend_str, str )==False): 
			print "ERROR: prepend_str=", prepend_str
			error_exit_with_message("prepend_str is not a string!")

		if(isinstance( append_str, str )==False): 
			print "ERROR: append_str=", append_str
			error_exit_with_message("append_str is not a string!")

		if(append_str=="" and prepend_str==""): error_exit_with_message('append_str=="" and prepend_str==""')


		sampler_post_process_filename = dirname(output_foldername) + '/' + prepend_str + basename(output_foldername).lower() + append_str
	
		return sampler_post_process_filename


#############################################################################################################

def get_sampler_post_process_filter_outfile(output_foldername):

		return get_sampler_post_process_generic_file(output_foldername , prepend_str="", append_str="_sample_filtered.out")


#############################################################################################################

def get_sampler_post_process_cat_outfile(output_foldername):

		return get_sampler_post_process_generic_file(output_foldername ,prepend_str="", append_str="_sample_catfile.out")

#############################################################################################################

