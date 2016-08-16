#!/usr/bin/env python

from sys import argv,exit
import sys
import traceback
from glob import glob
import string
import os
from os.path import basename, dirname, exists, expanduser
from time import sleep
import copy
######################################################################

from SWA_dagman_python.utility.SWA_util import *
from SWA_parse_options import parse_options, ensure_no_duplicate_options, get_option_name_args_safe

##############################################################################################################################
def benchmark_table_get_full_path(short_folder_name): #Can be a file as well!!!

	full_folder_name=""

	if(short_folder_name.count('Sept_2011_Ext_HD:')>0):
		full_folder_name=short_folder_name.replace('Sept_2011_Ext_HD:', '/Volumes/Sept_2011_Ext_HD_Parin/minirosetta/')

	elif(short_folder_name.count('Jan_11_Ext_HD:')>0):
		full_folder_name=short_folder_name.replace('Jan_11_Ext_HD:', '/Volumes/Jan_2011_Ext_HD_Parin/minirosetta/')

	elif(short_folder_name.count('Sept_10_Ext_HD:')>0):
		full_folder_name=short_folder_name.replace('Jan_11_Ext_HD:', '/Volumes/Parin_Ext_Sept_2010/minirosetta/')

	elif( short_folder_name[0:12]=='CS_DATA_SRC:' ):
		full_folder_name='/Users/sripakpa/minirosetta/Rosetta_rna_file/CS_PAPER/BENCHMARK_TABLE/DATA_SRC/' + short_folder_name[12:]

	elif( short_folder_name.count('/Users/sripakpa/minirosetta/')==0): 
		full_folder_name='/Users/sripakpa/minirosetta/' + short_folder_name

	else: 
		error_exit_with_message("Invalid short_foldername (%s) " %(shoft_folder_name))

	if(exists(full_folder_name)==False): error_exit_with_message("full_folder_name (%s) doesn't exist!" %(full_folder_name))

	#if(full_folder_name!=(os.path.abspath(full_folder_name)+'/')):
	#	print "full_folder_name                 = %s" %(full_folder_name)
	#	print "os.path.abspath(full_folder_name)= %s" %(os.path.abspath(full_folder_name)+'/')
	#	error_exit_with_message("full_folder_name!=(os.path.abspath(full_folder_name)+'/')")

	return full_folder_name


##############################################################################################################################
def get_benchmark_table_BMRB_file(job):

	if( job.has_key("src_folder")==False ): error_exit_with_message('job.has_key("src_folder")==False )')

	if( job.has_key("BMRB_file")==False ): error_exit_with_message('job.has_key("BMRB_file")==False )')

	BMRB_file=job["src_folder"] + job['BMRB_file']

	if(exists(BMRB_file)==False): error_exit_with_message("BMRB_file (%s) doesn't exist!" %(BMRB_file))

	return BMRB_file

##############################################################################################################################
def get_benchmark_table_native_pdb(job):

	if( job.has_key("src_folder")==False ): error_exit_with_message('job.has_key("src_folder")==False )')

	if( job.has_key("native_pdb")==False ): error_exit_with_message('job.has_key("native_pdb")==False )')


	#return job["native_pdb"] #Before Nov 05, 2011...used to return the native_pdb in the actual BIOX submit folder!

	native_pdb=job["src_folder"] + job["native_pdb"]

	if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!" %(native_pdb))

	return native_pdb

##############################################################################################################################

def get_benchmark_table_common_args_file_and_src_folder(job):

	common_args_file=""

	if(job.has_key("common_args") and job.has_key("src_folder") ): error_exit_with_message('job.has_key("common_args") and job.has_key("src_folder")')

	src_folder=""

	if(job.has_key("common_args") ): #Old (pre Nov 05, 2011) method!

		src_folder=job['common_args']

	elif( job.has_key("src_folder") ): #New (after Nov 05, 2011)!

		src_folder=job["src_folder"]
		
	else:
		error_exit_with_message("job does not have either the key 'common_args' or the key 'src_folder'")

	inner_common_args_name="COMMON_ARGS/common_args_region_FINAL.out"

	if(src_folder[-8:]=="$FARFAR$"):
		src_folder=src_folder[:-8]
		inner_common_args_name="FARFAR_common_args.txt"

	if(src_folder[-6:]=="$SELF$"):
		src_folder=src_folder[:-6]
	else:
		src_folder="%s/%s/" %(src_folder, job["folder_name"])

	src_folder=benchmark_table_get_full_path(src_folder)

	common_args_file="%s/%s" %(src_folder, inner_common_args_name)

	if(exists(src_folder)==False): error_exit_with_message("src_folder (%s) doesn't exist!" %(src_folder))

	if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!: " %(common_args_file)) 

	return (common_args_file, src_folder) 


#####################################################
def parse_motif_file(motif_file):

	if(exists(motif_file)==False): error_exit_with_message("motif_file (%s) doesn't exist!" %(motif_file))

	motif_list=[]
	motif_string_list=open( motif_file).readlines()
	for motif_string in motif_string_list:
		if(motif_string[-1]=="\n"): motif_string=motif_string[:-1]
		if(motif_string==""): continue
		if(motif_string[0]=="#"): continue
		motif_list.append(motif_string.split()[0]) # the split() get rid of extra space and tabs characters after the motif/folder name

	print '\n', "motif_list: " , motif_list, '\n'


	return motif_list

###############################################################
def parse_benchmark_data_source_file(data_src_file, star_mode, parse_silent_file): #This is meant to replace parse_benchmark_table_data_file()

	algorithm_list=[]
	data_src_list=[]

	line_list=open( data_src_file).readlines()

	data_src_string_list=[]

	abbreviation_map={}

	for line in line_list:

		line_split=line.split()

		if(len(line_split)==0): continue
		if(line_split[0][0]=="#"): continue

		if(line_split[0]=="ABBREVIATION"):
			
			if(len(line_split)!=2): error_exit_with_message("line_split[0]==\"COMMON_OPTIONS\" but len(line_split)!=2")

			if(line_split[1].count("=")!=1): error_exit_with_message("line_split[1].count(\"=\")!=1")	

			abbreviation_map["ABBREV:" +line_split[1].split("=")[0]]=line_split[1].split("=")[1]

		else:
			data_src_string_list.append(line_split)


	column_list=data_src_string_list[0]

	print '\n', column_list, '\n'


	if(column_list.count('folder_name')!=1): error_exit_with_message("column_list.count('folder_name')!=1") 
	if(column_list.count('main_folder')!=1): error_exit_with_message("column_list.count('main_folder')!=1")
	if(column_list.count('PRE_minimize_silent_file')!=1): error_exit_with_message("column_list.count('PRE_minimize_silent_file')!=1")
	if(column_list.count('silent_file_folder')!=1): error_exit_with_message("column_list.count('silent_file_folder')!=1")
	if(column_list.count('silent_file_type')!=1): error_exit_with_message("column_list.count('silent_file_type')!=1")
	if(column_list.count('old_SWA_idealize_helix')!=1): error_exit_with_message("column_list.count('old_SWA_idealize_helix')!=1")
	if(column_list.count('algorithm')!=1): error_exit_with_message("column_list.count('algorithm')!=1")
	if(column_list.count('common_args')!=1): error_exit_with_message("column_list.count('common_args')!=1")

	prev_data_src={}
	first_entry=True

	for n in range(1,len(data_src_string_list)): #exclude first column, which is the name column
		data_src_string=data_src_string_list[n]

		if(len(data_src_string)==0): continue #empty line
		if(data_src_string[0][0]=="#"): continue

		if(len(data_src_string)!=len(column_list)): 
			print "len(data_src_string) =%s | data_src_string=%s " %(len(data_src_string),  list_to_string(data_src_string))
			print "len(column_list)     =%s | column_list=%s " %( len(column_list), column_list)
			error_exit_with_message("len(data_src_string)!=len(column_list)") 

		data_src={} 

		for ii in range(len(column_list)): 
			column_name=column_list[ii]
			if(data_src_string[ii]=='""'):
				if(first_entry==True): error_exit_with_message('data_src_string[%d]="" but first_entry==True , data_src_string=%s ' %(ii, list_to_string(data_src_string) ) ) 
				if(column_name=="folder_name"): error_exit_with_message(' data_src_string[%d]=="" but column_name=="folder_name"!' %(ii)) 
				data_src[column_name]=prev_data_src[column_name]
			else:
				if(abbreviation_map.has_key(data_src_string[ii])):
					data_src_string[ii]=data_src_string[ii].replace(data_src_string[ii], abbreviation_map[data_src_string[ii]])

				data_src[column_name]=data_src_string[ii]

		if(data_src['folder_name'][0]=='#'): continue #redundant?

		first_entry=False

		if(parse_silent_file):
			if(data_src['silent_file_folder']=="MISSING"): continue  #Important to have this statement before update prev_data_src
			if(data_src['silent_file_type']=="MISSING"): continue    #Important to have this statement before update prev_data_src

		prev_data_src=copy.deepcopy(data_src)

		if(star_mode==True): 
			if( data_src['folder_name'][0]!='*'): continue
	
		if( data_src['folder_name'][0]=='*'): data_src['folder_name']=data_src['folder_name'][1:]

		if(data_src['algorithm'] not in algorithm_list): algorithm_list.append(data_src['algorithm'])

		data_src['main_folder']=benchmark_table_get_full_path(data_src['main_folder'])

		if(parse_silent_file):
			data_src['silent_file_folder']=benchmark_table_get_full_path(data_src['silent_file_folder'])
		else:
			data_src['silent_file_folder']=""
			data_src['silent_file_type']=""

		print "data_src #%d " %(n), data_src

		data_src_list.append(data_src)

	algorithm_list.sort()

	return (data_src_list, algorithm_list)

###############################################################
def parse_benchmark_table_data_file(job_list_file, algorithm, star_mode):

	algorithm_list=[]
	job_list=[]

	job_string_list=open( job_list_file).readlines()

	for n in range(len(job_string_list)):
		job_string_list[n]=job_string_list[n].split()

	column_list=job_string_list[0]

	print column_list

	if(column_list.count('folder_name')!=1): error_exit_with_message("column_list.count('folder_name')!=1") 
	if(column_list.count('main_folder')!=1): error_exit_with_message("column_list.count('main_folder')!=1")
	if(column_list.count('silent_file')!=1): error_exit_with_message("column_list.count('silent_file')!=1")
	if(column_list.count('rmsd_col_name')!=1): error_exit_with_message("column_list.count('rmsd_col_name')!=1")
	if(column_list.count('algorithm')!=1): error_exit_with_message("column_list.count('algorithm')!=1")
	if(column_list.count('PRE_silent_file')!=1): error_exit_with_message("column_list.count('PRE_silent_file')!=1")



	if(column_list.count('common_args')!=1 and column_list.count('src_folder')!=1): 
		error_exit_with_message("column_list.count('common_args')!=1 and column_list.count('src_folder')!=1")

	if(column_list.count('native_pdb')!=1):  error_exit_with_message("column_list.count('native_pdb')!=1")


	prev_job={}
	first_job=True

	for n in range(1,len(job_string_list)): #exclude first column, which is the name column
		job_string=job_string_list[n]

		if(len(job_string)==0): continue #empty line
		if(job_string[0][0]=="#"): continue

		if(len(job_string)!=len(column_list)): 
			print "job_string= ", job_string
			print "column_list= ", column_list
			error_exit_with_message("len(job_string)!=len(column_list)") 

		job={} 

		for ii in range(len(column_list)): 
			column_name=column_list[ii]
			if(job_string[ii]=='""'):
				if(first_job==True): error_exit_with_message('job_string[%d]="" but first_job==True , job_string=%s ' %(ii, list_to_string(job_string) ) ) 
				if(column_name=="folder_name"): error_exit_with_message(' job_string[%d]=="" but column_name=="folder_name"!' %(ii)) 
				job[column_name]=prev_job[column_name]
			else:
				job[column_name]=job_string[ii]

		if(job['folder_name'][0]=='#'): continue

		first_job=False

		prev_job=copy.deepcopy(job)

		if(star_mode==True): 
			if( job['folder_name'][0]!='*'): continue
	
		if( job['folder_name'][0]=='*'): job['folder_name']=job['folder_name'][1:]

		if(job['algorithm'] not in algorithm_list): algorithm_list.append(job['algorithm'])

		if(algorithm=="create_table"):
			if(job['silent_file']=="MISSING"): continue
			job['PRE_silent_file']="BLAH_PRE_silent_file_BLAH"
			job['PRE_silent_file_list']="BLAH_PRE_silent_file_list_BLAH"



		if(algorithm=="recal_rmsd"): 
			if(job['PRE_silent_file']=="MISSING"): continue
			#if(job['common_args']=="MISSING"): continue Comment out this line on Nov 06, 2011
			job['silent_file']=="BLAH_silent_file_BLAH"
			job['silent_file_list']="BLAH_silent_file_list_BLAH"


		job['main_folder']=benchmark_table_get_full_path(job['main_folder'])

		(job['common_args'], job['src_folder'])=get_benchmark_table_common_args_file_and_src_folder(job)

		job['native_pdb']=get_benchmark_table_native_pdb(job)

		if( job.has_key('BMRB_file') ): job['BMRB_file']=get_benchmark_table_BMRB_file(job)

		print job, "job #%d " %(n)

		job_list.append(job)

	algorithm_list.sort()

	return (job_list, algorithm_list)

############################################################################################

def parse_benchmark_job_file(job_list_file, setup_folders_and_files, star_mode, verbose):

	HOMEDIR = expanduser('~') 

	PRISTINE_FOLDER=HOMEDIR + "/minirosetta/Rosetta_rna_file/SWA_LOOP_PAPER/"

	if(exists(job_list_file)==False): error_exit_with_message("job_list_file (%s) doesn't exist!!" %(job_list_file))  

	line_list=open( job_list_file).readlines()
	job_string_list=[]
	common_option_list=[]

	for line in line_list:

		line_split=line.split()

		if(len(line_split)==0): continue
		if(line_split[0][0]=="#"): continue

		if(line_split[0]=="COMMON_OPTIONS"):
			
			if(len(line_split)<=4): error_exit_with_message("line_split[0]==\"COMMON_OPTIONS\" but len(line_split)<=4")
			if(line_split[1]!='"' or line_split[-1]!='"') : error_exit_with_message('line_split[0]==\"COMMON_OPTIONS\" but line_split[2]!=="" or job_string[-1]!==""')
			common_option_list=line_split[2:-1]
			if(verbose): print "common_option_list=%s" %(list_to_string(common_option_list))
			if(len(common_option_list)==0): error_exit_with_message("len(common_option_list)==0!")

		else:
			job_string_list.append(line_split)


	column_list=job_string_list[0]

	if(verbose): print column_list

	if(column_list.count('folder_name')!=1): error_exit_with_message("column_list.count('folder_name')!=1") 
	if(column_list.count('native_pdb')!=1): error_exit_with_message("column_list.count('native_pdb')!=1") 
	if(column_list.count('cutpoint_open')!=1): error_exit_with_message("column_list.count('cutpoint_open')!=1") 
	if(column_list.count('bulge_res')!=1): error_exit_with_message("column_list.count('bulge_res')!=1") 

	job_list=[]

	total_slave_nodes=0

	for n in range(1,len(job_string_list)): #exclude first column, which is the name column

		job_string=job_string_list[n]

	#	if(len(job_string)!=len(column_list)): error_exit_with_message("len(job_string)!=len(column_list)") 

		job={} 

		for ii in range(len(column_list)): 
			column_name=column_list[ii]
			if(ii!=(len(column_list)-1)):

				if(column_name=="FARFAR_nstruct"): error_exit_with_message("column_name==\"FARFAR_nstruct\". Please instead specify nstruct as a direct option to create_benchmakr_job_files.py")
				if(column_name=="SWA_nstruct"): 	 error_exit_with_message("column_name==\"SWA_nstruct\". Please instead specify nstruct as a direct option to create_benchmakr_job_files.py")
				if(column_name=="nstruct"):			   error_exit_with_message("column_name==\"nstruct\". Please instead specify nstruct as a direct option to create_benchmakr_job_files.py")


				job[column_name]=job_string[ii]

			else: #the extra column
				if(column_name!="zextra"):  error_exit_with_message('last_column_name!="extra"') 
				if(len(job_string)<(ii+2) ): error_exit_with_message('len(job_string)<(ii+2)"')
				if(job_string[ii]!='"' or job_string[-1]!='"') : error_exit_with_message('job_string[ii]!=="" or job_string[-1]!==""')
	
				if(len(job_string)==(ii+2)):
					job["zextra"]=""			
				else:
					job["zextra"]=list_to_string(job_string[ii+1:-1])

				############################Oct 18, 2011, Add common_options#############
				job["zextra"]=job["zextra"] + " " + list_to_string(common_option_list)

				ensure_no_duplicate_options(job["zextra"], column_list)

		#if(job['folder_name'][0]=='#'): continue , this is now taken care of above...

		PRISTINE_motif_folder_name=""

		if(star_mode==True): 
			if( job['folder_name'][0]!='*'): continue
	
		if( job['folder_name'][0]=='*'): job['folder_name']=job['folder_name'][1:]

		if(dirname(job['folder_name'])==""):
			PRISTINE_motif_folder_name=PRISTINE_FOLDER + job['folder_name']
		else:
			if(job['folder_name'][0]!="/"): error_exit_with_message('dirname(job[\'folder_name\'])!="" but job[\'folder_name\'][0]!="/", job[\'folder_name\']=%s' %(job['folder_name']) ) 
			PRISTINE_motif_folder_name=HOMEDIR + "/minirosetta/Rosetta_rna_file/" + job['folder_name'] #os.path.abspath()
			job['folder_name']=basename(job['folder_name'])

		if(exists(PRISTINE_motif_folder_name)==False): error_exit_with_message("PRISTINE_motif_folder_name (%s) doesn't exist!" %(PRISTINE_motif_folder_name))

		job['PRISTINE_motif_folder_name']=PRISTINE_motif_folder_name

		bulge_string_list=job['bulge_res'].split('/')

		if(len(bulge_string_list)>2): error_exit_with_message("len(bulge_string_list(%s))>2" %(list_to_string(bulge_string_list) ) ) 

		for n in range(len(bulge_string_list)):

			bulge_res_string=bulge_string_list[n]

			#Shallow copies of dictionaries can be made using dict.copy(), and of lists by assigning a slice of the entire list, for example, copied_list = original_list[:].

			job_act=copy.deepcopy(job)

			job_act['folder_name']=job['folder_name']
			job_act['bulge_res']=""
			split_string=bulge_res_string.split("-")

			for n in range(len(split_string)):
				seq_num=int(split_string[n]) # the purpose of this is check that the format is valid
				if(n!=0): job_act['bulge_res']+="-"
				job_act['bulge_res']+="%d" %(seq_num)			

			if(setup_folders_and_files):
				submit_subprocess("mkdir %s " %job_act['folder_name'])			

				native_pdb_location=PRISTINE_motif_folder_name + "/" + job_act['native_pdb']
				if(exists(native_pdb_location)==False): error_exit_with_message("native_pdb %s doesn't exist!" %(native_pdb_location) )
				submit_subprocess("cp %s %s/%s " %(native_pdb_location	, job_act['folder_name'] , job_act['native_pdb']	 ) )

				extra_pdb_list=job_act['zextra'].split()

				for extra_pdb in extra_pdb_list:
					if((extra_pdb.count(".pdb")!=1) and (extra_pdb.count(".bmrb")!=1)): continue #check if it a pdb file

					pdb_location=PRISTINE_motif_folder_name + "/" + extra_pdb
					if(exists(pdb_location)==False): error_exit_with_message("pdb %s doesn't exist!" %(pdb_location) ) 
					submit_subprocess("cp %s %s/%s " %(pdb_location	, job_act['folder_name'],  extra_pdb ) )
	
			if(verbose): print " job #", (len(job_list)+1), " ", job_act

#Switch to above since below doesn't work in python 2.3 (Biox2)
#			print " job #%d {" %(len(job_list)+1),
#			for key in sorted(job_act.iterkeys()):
#				print "'%s': '%s'" % (key, job_act[key]),
#			print "}"

			total_slave_nodes+=int(job_act['num_slave_nodes'])

			job_list.append(job_act)

	if(verbose): 
		print "------------------------------------------------------------------------------------------"
		print "total_slave_nodes REQUIRED TO RUN ALL JOBS: %d" %(total_slave_nodes) 
		print "------------------------------------------------------------------------------------------"

	return job_list


#####################Functions below could below to a new SWA_parse_benchmark_util.py###################

def find_matching_job(job_submission_list,motif):

	found_job=False
	job={}

	for curr_job in job_submission_list:
		if(curr_job['folder_name']==motif): 
			if(found_job): error_exit_with_message("Already found_job for motif (%s)" %(motif))
			job=copy.deepcopy(curr_job)
			found_job=True

	return (job, found_job)

##########################################################################
def create_fasta_from_data_src(data_src, native_pdb):

	common_args_file= "%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], data_src['common_args'])

	if(exists(native_pdb)==False): error_exit_with_message("native_pdb (%s) doesn't exist!: " %(native_pdb) ) 
	if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file))

	common_args_string=safe_readlines( common_args_file )[0]
	fasta_file=get_option_name_args_safe(common_args_string, "fasta", "")
	submit_subprocess("SWA_pdb2fasta.py %s > %s" %(native_pdb, fasta_file))


##########################################################################
def create_VDW_rep_pdb_from_data_src(data_src):

	common_args_file= "%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], data_src['common_args'])

	if(exists(common_args_file)==False): error_exit_with_message("common_args_file (%s) doesn't exist!" %(common_args_file))

	common_args_string=safe_readlines( common_args_file )[0]

	VDW_rep_screen_info=get_option_name_args_safe(common_args_string, "VDW_rep_screen_info", [""])

	if(len(VDW_rep_screen_info) % 3 != 0): error_exit_with_message("len(VDW_rep_screen_info) % 3 != 0")

	for n in range(len(VDW_rep_screen_info)):

		if(n % 3 != 0): continue

		VDW_rep_screen_PDB="%s/%s/%s" %(data_src['main_folder'], data_src['folder_name'], VDW_rep_screen_info[n])

		if(exists(VDW_rep_screen_PDB)==False): error_exit_with_message("VDW_rep_screen_PDB (%s) doesn't exist!" %(VDW_rep_screen_PDB))

		submit_subprocess("cp %s %s" %(VDW_rep_screen_PDB, basename(VDW_rep_screen_PDB) ) )

##########################################################################

