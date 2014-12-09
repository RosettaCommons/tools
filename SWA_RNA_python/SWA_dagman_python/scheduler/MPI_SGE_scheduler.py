#!/usr/bin/env python

######################################################################
###!/usr/bin/env python2.7-mpi ##python2.7-mpi should be invoke in DAG_continuous_MPI_wrapper.py, not here!
###OK doesn't make sense to have /usr/bin/env python2.7-mpi here, since this script is never called as an exeutable!!

#DO NOT INCLUDE 'from mpi4py import MPI' line here! Since on some machine such as LONESTAR TACC. This will lead to error if standard python is used and not python2.7-mpi.
#PROBLEM IS THAT THIS SCRIPT IS CALLED BY DAG_continuous_MPI_wrapper.py which uses python2.7-mpi and submit_DAG_job.py which uses normal python.

from SWA_dagman_python.utility.SWA_util import *

from scheduler_common import *

MPI_SCHEDULER=True

SCHEDULER_TYPE="MPI_SGE_scheduler.py"

######################################################################


#Use this to run cluster using a SGE (Sun Grid Engine) scheduler in MPI mode. 

#Shown below is a table of common job script options in PBS (TORQUE) and their SGE equivalents. Directives in SGE are marked with a #$ rather than the #PBS used by PBS.
# PBS 	SGE
# -l walltime=hours:minutes:seconds 	 -l h_rt=<hours:minutes:seconds>
# -l pmem=<size> 	 -l mem=<size>
# -N <name> 	 -N <name>
# -A <project> 	 -P <project>
# -l nodes=<n>:ppn=<np> 	 -pe <environment> <number> 

#     qsh    -  submit an interactive  X-windows  session  to  Sun
#               Grid Engine.


#MPI assumes, implicitly, the existence of an environment in which an application runs. It does not provide ``operating system'' services, such as a general ability to query what processes are running, to kill arbitrary processes, to find out properties of the runtime environment (how many processors, how much memory, etc.). 

######################################################################

def running_job_state_str():
	return 'R'

def pending_job_state_str():
	return 'Q' #job is queued, eligible to run or routed. 


######################################################################

def update_scheduler_queue_log_command(): #For MPI, currently do nothing.

	pass

######################################################################

def get_queued_jobs_status():

	from mpi4py import MPI ##NEED THIS!
	#This should be ok. A module can contain executable statements as well as function definitions. These statements are intended to initialize the module. They are executed only the first time the module is imported somewhere. [1] (http://docs.python.org/tutorial/modules.html#id2)

	MPI_comm = MPI.COMM_WORLD
	MPI_rank = MPI_comm.Get_rank() #Numbering of current core
	MPI_size = MPI_comm.Get_size() #Total number of cores

	MPI_name = MPI.Get_processor_name()

	if(MPI_rank!=0): error_exit_with_message("Only the master_process should call this function but MPI_rank!=0 !!!")

	job_info_list=[]

	for slave_job_num in range(0, MPI_size-1): #range is 0 to MPI_size-2 (or equivalently num_slave_jobs-1)

		job_info={}

		job_info['JOBID']=slave_job_num

		job_info['JOB_NAME']="%s_%s" %(abspath( 'SLAVE_JOBS/' ).replace('/','_')[1:], slave_job_num)  #Assume specific form!

		job_info['USER']="BLAH"
		
		job_info['EXEC_HOST']=MPI_name 

		job_info['STATE']=running_job_state_str()

		job_info_list.append(job_info)

	'''
	job_info['JOBID'], job_info['JOB_NAME'], job_info['USER'], job_info['EXEC_HOST'], job_info['STATE']
	'''

	return job_info_list



######################################################################
def 	master_kill_all_slave_jobs_and_exit_scheduler_specific(exit_message):

	from mpi4py import MPI ##NEED THIS!
	#This should be ok. A module can contain executable statements as well as function definitions. These statements are intended to initialize the module. They are executed only the first time the module is imported somewhere. [1] (http://docs.python.org/tutorial/modules.html#id2)

	print "A error had occured. About to abort MPI job, error_message= " + exit_message

	'''
	#Comment out this section. Problem is not the compute node might not have permission to issue qdel command (this is true for sure on Lonestar)
	MPI_JOB_ID_filename='MPI_JOB_ID.txt'

	if(exists(MPI_JOB_ID_filename)==False): error_exit_with_message("MPI_JOB_ID_filename (%s) doesn't exist!" %(MPI_JOB_ID_filename))

	MPI_JOB_ID_lines=safe_open(MPI_JOB_ID_filename, 'r').readlines() 

	if(len(MPI_JOB_ID_lines)!=1): error_exit_with_message("len(MPI_JOB_ID_lines)=(%s)!=1" %(len(MPI_JOB_ID_lines)))

	MPI_JOB_ID_first_line_list=MPI_JOB_ID_lines[0].split()

	if(len(MPI_JOB_ID_first_line_list)!=1): error_exit_with_message("len(MPI_JOB_ID_lines[0].split())=(%s)!=1" %(MPI_JOB_ID_first_line_list) )

	MPI_JOB_ID=MPI_JOB_ID_first_line_list[0]
	
	print "MPI_JOB_ID= (%s)" %(MPI_JOB_ID)	

	MPI_JOB_ID=int(MPI_JOB_ID)

	command='qdel %d ' %(MPI_JOB_ID)

	print command
	submit_subprocess( command )
	'''

	MPI.COMM_WORLD.Abort() 
	
	#If use qsub, then once MPI is aborted, then all child subprocess (e.g. Rosetta execution should dies right away). Note this is not true if running MPI on local machine or server's frontend 
	#Note I confirmed this for native_screen on the Lonestar by checking that the time-stamp of all Rosetta execution log_file is prior or same time as the time-stamp of MPI_log.out

######################################################################
def queue_job_command(job_name, outfile, errfile, job_script, job_dir_name, walltime, memory_limit, memory_reserve, queue_name):

	#In the MPI paradigm, all process are submitted together not individually

	error_exit_with_message("queue_job_command() function not implemented for MPI scheduler")



######################################################################
'''
login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ qstat -r
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
 602986 0.56676 work_01462 parins       qw    06/06/2012 11:09:55                                   48        
       Full jobname:     work_01462_parins_June_06_2012_SWA_GAAA_tetraloop_lonestar_test_TRAIL_3_LONESTAR_XSEDE_TEST_MPI
       Requested PE:     12way 48
       Hard Resources:   mem_free=2048M (0.000000)
                         mem_total=23.4G (0.000000)
                         h_rt=86400 (0.000000)
       Soft Resources:   
       Hard requested queues: normal
 602991 0.00000 work_01462 parins       qw    06/06/2012 11:12:09                                   48        
       Full jobname:     work_01462_parins_June_06_2012_SWA_GAAA_tetraloop_lonestar_test_TRAIL_3_LONESTAR_XSEDE_TEST_MPI
       Requested PE:     12way 48
       Hard Resources:   mem_free=2048M (0.000000)
                         mem_total=23.4G (0.000000)
                         h_rt=86400 (0.000000)
       Soft Resources:   
       Hard requested queues: normal

login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
 602986 0.56789 work_01462 parins       qw    06/06/2012 11:09:55                                   48        
 602991 0.53418 work_01462 parins       qw    06/06/2012 11:12:09                                   48        
login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ qdel 602986
parins has deleted job 602986
login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ qstat
job-ID  prior   name       user         state submit/start at     queue                          slots ja-task-ID 
-----------------------------------------------------------------------------------------------------------------
 602991 0.53367 work_01462 parins       qw    06/06/2012 11:12:09                                   48        
login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ 

login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ qdel -u parins
parins has deleted job 602991
login2:/work/01462/parins/June_06_2012_SWA_GAAA_tetraloop_lonestar_test/TRAIL_3_LONESTAR_XSEDE_TEST$ qstat
'''

######################################################################
def submit_DAG_job_scheduler_specific(master_wall_time, master_memory_reserve, num_slave_nodes, dagman_file, verbose):

	wall_time=master_wall_time #In MPI implementation, no distinction between master_wall_time and slave_wall_time

	memory_reserve=master_memory_reserve #In MPI implementation, no distinction between master_memory_reserve and memory_reserve

	#walltime (is hours): the job wall-clock run time limit.

	#memory_reserve (in MB): Amount of memory reserved/allocated for this job. This parameter guarantee memory allocation [Note memory_reserve should be less-than-or-equal-to memory_limit].

	#Always set memory_limit to equal memory_reserve.
	#memory_limit (in MB): job is killed when it exceeds this memory limit. This parameter does not guarantee memory allocation, it's just a threshold. This parameter is not implemented for some scheduler such as Portable Batch System (PBS).

	#Assume using default queue. (For BIOX2 cluster, this is the SP qeuue | #PBS -q SP)

	MPI_DAG_script_EXE=get_PYEXE("dagman/DAG_continuous_MPI_wrapper.py")

	MPI_job_name = abspath( "MPI" ).replace("/","_")[1:]

	MPI_outfile='MPI_log.out'
	MPI_errfile='MPI_log.err'

	if( exists(MPI_outfile) ): submit_subprocess("rm %s" %(MPI_outfile))
	if( exists(MPI_errfile) ): submit_subprocess("rm %s" %(MPI_errfile))
	#NOTE: the master_logs are taken care off by DAG_continuous_MPI_wrapper.py

	if( exists("SLAVE_JOBS/") ): submit_subprocess("rm -rf SLAVE_JOBS/")

	#Don't actually need the -j option, since num_slave_nodes=MPI_comm.Get_size()-1.
	MPI_job_script=" %s -num_slave_jobs %d -dagman_file %s" %(MPI_DAG_script_EXE, num_slave_nodes, dagman_file)

	COMPUTER_CLUSTER_NAME=str(os.environ['COMPUTER_CLUSTER_NAME'])

	CORES_PER_NODE=0

	MPI_QUEUE_NAME="BLAH"

	if(COMPUTER_CLUSTER_NAME=="LOCAL_TSET"):

		CORES_PER_NODE=1
		MPI_QUEUE_NAME="normal"

	elif((COMPUTER_CLUSTER_NAME=="LONESTAR-TACC-XSEDE")):

		CORES_PER_NODE=12 #This is specific to the architect of the LONESTAR cluster. Each node contains 12 cores and each submission must use the entire node(s).
		MPI_QUEUE_NAME="normal" 	#"development" #This is specific to the Lonestar cluster

		#if(num_slave_nodes <=23 and wall_time<=1): MPI_QUEUE_NAME="development"
			

		#Get a list of the queues you can submit stuff to: qconf -sql  (Sun Grid Engine)
		'''LONSTAR XSEDE.
		login1:/work/01462/parins/minirosetta/TRAIL_2_June_06_2012_SWA_GAAA_tetraloop_lonestar_test/AQUA_1_LONESTAR_TEST$ qstat -g c
		CLUSTER QUEUE                   CQLOAD   USED    RES  AVAIL  TOTAL aoACDS  cdsuE  
		--------------------------------------------------------------------------------
		development                       0.72    636      0     12    768    756      0 
		gpu                               0.13     72      0     24     96     72      0 
		grace                             0.00      0      0     72     72      0      0 
		grace-serial                      0.00      0      0     12     12      0      0 
		largemem                          0.00      0      0    336    336      0      0 
		normal                            0.80  20832      0    480  22464  20724     12 
		request                           0.80      0      0   1728  22464  20724     12 
		serial                            0.73    192      0      0    768    768      0 
		stci                              0.70      0      0      0    192    192      0 
		sysdebug                          0.80      0      0   1752  23232  21480      0 
		systest                           0.80      0      0   1752  23232  21480      0 
		vis                               0.13      0      0     12     96     72     12 
		login1:/work/01462/parins/minirosetta/TRAIL_2_June_06_2012_SWA_GAAA_tetraloop_lonestar_test/AQUA_1_LONESTAR_TEST$ 
		'''


	elif((COMPUTER_CLUSTER_NAME=="RANGER-TACC-XSEDE")):

		CORES_PER_NODE=16 #This is specific to the architect of the RANGER cluster. Each node contains 16 cores and each submission must use the entire node(s).
		MPI_QUEUE_NAME="normal" 	#"development" #This is specific to the Lonestar cluster

	elif COMPUTER_CLUSTER_NAME == "STAMPEDE-TACC-XSEDE":

		# This is specific to the architect of the STAMPETE cluster.
		CORES_PER_NODE = 16

		MPI_QUEUE_NAME ="normal"

	else:
		error_exit_with_message("Unsupported COMPUTER_CLUSTER_NAME=%s" %(COMPUTER_CLUSTER_NAME))

	print "ENVIROMENT_VARIABLE: COMPUTER_CLUSTER_NAME=%s | CORES_PER_NODE=%s | MPI_QUEUE_NAME=%s " %(COMPUTER_CLUSTER_NAME, CORES_PER_NODE, MPI_QUEUE_NAME)

	#E.G: Requests X (=12) tasks/node, Y (=24) cores total (Y must be multiples of 12, set X to 12 for lonestar)

	if( ( (num_slave_nodes+1) % CORES_PER_NODE ) != 0 ): 
		error_exit_with_message("num_slave_nodes+1 (%s) is not divisible by CORES_PER_NODE (%s)" %(num_slave_nodes+1, CORES_PER_NODE))
	

	#FROM: http://www.biostat.jhsph.edu/bit/cluster-usage.html
	#qsub -cwd -l mem_free=MEM_NEEDED,h_vmem=MEM_MAX  batch.sh
	#Note http://th.physik.uni-frankfurt.de/wiki-it/index.php?title=Sun_Grid_Engine used virtual_free in in place of mem_free.

	##$ -M yourmail@stanford.edu     # Address for email notification
	##$ -m be        # Email at Begin and End of job
	##	set -x  # Echo commands, use set echo with csh

	if COMPUTER_CLUSTER_NAME == "STAMPEDE-TACC-XSEDE":
		# Use sbatch instead of the standard qsub

		#Submit the MPI job. This includes both the master_job and all the slave_jobs.
		submit_file = "SBATCH_MPI_DAG_JOB.sh" 

		if exists(submit_file): 
			submit_subprocess("rm %s" % submit_file)

		JOB = open(submit_file, 'w')

		JOB.write('#!/bin/bash\n\n')

		#QSUB_JOB.write( '#$ -cwd\n\n' ) ## Start job in submission directory

		JOB.write('#SBATCH -J %s\n\n' % MPI_job_name) # Job Name
		JOB.write('#SBATCH -o %s\n\n' % (MPI_outfile + "_SBATCH")) # Out log_file
		JOB.write('#SBATCH -e %s\n\n' % (MPI_errfile + "_SBATCH")) # Err log_file
		JOB.write('#SBATCH -t %d:00:00\n\n' % wall_time)  # Run time (hh:mm:ss)

		# Request num_slave_nodes + 1 number of task. The job will automatically
		# acquices enough nodes to executes (num_slave_nodes + 1) tasks, with
		# On Stampede, the default is 16 tasks per node.
		# e.g. #SBATCH -n 64 would run on 4 nodes, 16 tasks per node.

		JOB.write('#SBATCH -n %s\n\n' % (num_slave_nodes + 1)) # number of tasks
		JOB.write('#SBATCH -N %s\n\n' % ((num_slave_nodes + 1) / CORES_PER_NODE) ) # number of nodes, redundunt.

		JOB.write('#SBATCH --mail-user=sripakpa@stanford.edu\n\n') # Address for email notification
		JOB.write('#SBATCH --mail-type=ALL\n\n')  # "ALL" will notify user of any state changes (BEGIN, END, FAIL, REQUEUE)

		# JOB.write( '#SBATCH -A TG-MCB090153\n\n') # OLD Das Lab Project Charge No.

		JOB.write('#SBATCH -A TG-MCB120152\n\n') # New Das Lab Project Charge No.

		# I couldn't get memory allocation option to work on STAMPEDE.
		'''
		if memory_reserve != 0:
			# minimum REAL memory required per CPU process in MegaBytes.
			JOB.write('#SBATCH --mem-per-cpu=%dM\n\n' % (memory_reserve))
		'''

		JOB.write('#SBATCH -p %s\n\n' % MPI_QUEUE_NAME) # queue to submit the job.

		JOB.write('#SBATCH --verbose\n\n') # Verbose!
		JOB.write('#SBATCH --verbose\n\n') # Further Increase Verbosity!

		JOB.write( 'echo $SLURM_JOB_ID > MPI_JOB_ID.txt\n\n')

		# ibrun is specific command for submitting MPI jobs TACC clusters.
		JOB.write( 'ibrun %s >%s 2>%s\n\n' %(MPI_job_script, MPI_outfile, MPI_errfile)) 

		JOB.close()

		submit_subprocess('sbatch %s' % submit_file)

	else:

		#Submit the MPI job. This includes both the master_job and all the slave_jobs.
		qsub_submit_file="QSUB_MPI_DAG_JOB.sh" 

		if(exists(qsub_submit_file)): submit_subprocess("rm %s" %(qsub_submit_file))

		QSUB_JOB = open( qsub_submit_file, 'w')

		QSUB_JOB.write( '#!/bin/bash\n\n' )

		QSUB_JOB.write( '#$ -cwd\n\n' ) ## Start job in submission directory

		QSUB_JOB.write( '#$ -N %s\n\n' %(MPI_job_name))
		QSUB_JOB.write( '#$ -o %s\n\n' %(MPI_outfile + "_QSUB" ))
		QSUB_JOB.write( '#$ -e %s\n\n' %(MPI_errfile + "_QSUB" ))
		QSUB_JOB.write( '#$ -l h_rt=%d:00:00\n\n' %(wall_time))  # Run time (hh:mm:ss)
		QSUB_JOB.write( '#$ -pe %dway %d \n\n' %(CORES_PER_NODE, num_slave_nodes+1) )   #e.g. -pe 12way 24 --> Requests 12 cores/node, 24 cores total

		QSUB_JOB.write( '#$ -M sripakpa@stanford.edu\n\n') # Address for email notification
		QSUB_JOB.write( '#$ -m bea\n\n')  # "bea" mean mail is sent at the beginning, end, and at abort time (if it happens) 

		QSUB_JOB.write( '#$ -A TG-MCB090153\n\n') # Das Lab Project Charge No.

		if(memory_reserve!=0):
			QSUB_JOB.write( '#$ -l mem_free=%dM\n\n' %(memory_reserve)) #memory per CPU process. 
			QSUB_JOB.write( '#$ -l h_vmem=%dM\n\n' %(memory_reserve)) #memory per CPU process. 


		QSUB_JOB.write( '#$ -q %s\n\n' %(MPI_QUEUE_NAME))

		QSUB_JOB.write( 'echo $JOB_ID > MPI_JOB_ID.txt\n\n')

		# ibrun is specific command for submitting MPI jobs TACC clusters.
		QSUB_JOB.write( 'ibrun %s >%s 2>%s\n\n' %(MPI_job_script, MPI_outfile, MPI_errfile))

		QSUB_JOB.close()

		submit_subprocess('qsub -V %s' %(qsub_submit_file)) #-V --> inherit the submission environment


