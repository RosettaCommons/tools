#!/usr/bin/env python

######################################################################
###!/usr/bin/env python2.7-mpi ##python2.7-mpi should be invoke in DAG_continuous_MPI_wrapper.py, not here!
###OK doesn't make sense to have /usr/bin/env python2.7-mpi here, since this script is never called as an exeutable!!

#DO NOT INCLUDE 'from mpi4py import MPI' line here! Since on some machine such as LONESTAR TACC. This will lead to error if standard python is used and not python2.7-mpi.
#PROBLEM IS THAT THIS SCRIPT IS CALLED BY DAG_continuous_MPI_wrapper.py which uses python2.7-mpi and submit_DAG_job.py which uses normal python.

from SWA_dagman_python.utility.SWA_util import *

from scheduler_common import *

MPI_SCHEDULER=True

SCHEDULER_TYPE="MPI_PBS_scheduler.py"

######################################################################

#Use this to run cluster using a PBS (TORQUE)) scheduler in MPI mode.. 

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

	NUM_NODES=0

	MPI_QUEUE_NAME="BLAH"

	PBS_REQUEST_COMPUTE_NODES_LINE=""
	MPI_EXEC_WRAPPER_STR=""

	if(COMPUTER_CLUSTER_NAME=="LOCAL_TSET"):

		CORES_PER_NODE=1
		NUM_NODES=1
		MPI_QUEUE_NAME="normal"

	elif((COMPUTER_CLUSTER_NAME=="TRESTLES-SDSC-XSEDE")):

		CORES_PER_NODE=32 #This is specific to the architect of the LONESTAR cluster. Each node contains 12 cores and each submission must use the entire node(s).
		MPI_QUEUE_NAME="normal" 	#"development" #This is specific to the Lonestar cluster


		if( ( (num_slave_nodes+1) % CORES_PER_NODE ) != 0 ): 
			error_exit_with_message("num_slave_nodes+1 (%s) is not divisible by CORES_PER_NODE (%s)" %(num_slave_nodes+1, CORES_PER_NODE))

		NUM_NODES=int((num_slave_nodes+1) / CORES_PER_NODE )

		#NOTE: PBS -l pmem=2048mb appears to made job stuck on QUEUE for TRESTLES. This is maybe becuase the memory capacity of TRESTLE node (32 core) is 64 GB

		#MAX_MEMORY_PER_CORE=1956   #95% of 2048M   (2048 * 0.95 = 1 945.6) #Comment this out on June 17, 2012. OK this works if NUM CORES is small (64 cores/ nodes). But for large number of cores (512 cores /16 nodes), the MPI run just freezes after start running.

		MAX_MEMORY_PER_CORE=0  #Switch to this on June 17, 2012

		PBS_REQUEST_COMPUTE_NODES_LINE='-l nodes=%d:ppn=%d' %(NUM_NODES, CORES_PER_NODE) 

		MPI_EXEC_WRAPPER_STR="mpirun -np %s -hostfile $PBS_NODEFILE " %(num_slave_nodes+1)


		'''Trestles XSEDE.
		Queue Name 	Limits 	Comments
		normal 	max wallclock = 48 hours
		max node count = 32 	Used for exclusive access to entire compute nodes (32 processors). This queue should be used for jobs that require at least 32 processors or large memory. Note: the max wallclock can be extended to up to 2 weeks. Send email to help@xsede.org to request this for your job.
		shared 	max wallclock = 48 hours
		max node count = 4 	Used for shared access to a node. This queue is useful for debugging codes or for jobs that require less than 32 processors/node. Note: the max wallclock can be extended to up to 2 weeks. Send email to help@xsede.org to request this for your job.


		trestles-login2:~/minirosetta/TRAIL_1_June_10_2012_SWA_GAAA_tetraloop_lonestar_test/AQUA_2_TRESTLE_SHARE_QUEUE$ qstat -q

		server: trestles-fe1

		Queue            Memory CPU Time Walltime Node  Run Que Lm  State
		---------------- ------ -------- -------- ----  --- --- --  -----
		shared             --      --       --      --   63   2 --   E R
		normal             --      --       --      --   26   8 --   E S
				                                           ----- -----
				                                              89    10
		trestles-login2:~/minirosetta/TRAIL_1_June_10_2012_SWA_GAAA_tetraloop_lonestar_test/AQUA_2_TRESTLE_SHARE_QUEUE$ 
		'''


	elif((COMPUTER_CLUSTER_NAME=="KRAKEN-NICS-XSEDE")):

		CORES_PER_NODE=12
		MPI_QUEUE_NAME=""

		MAX_MEMORY_PER_CORE=0

		PBS_REQUEST_COMPUTE_NODES_LINE='-l size=%s' %(num_slave_nodes+1) 

		MPI_EXEC_WRAPPER_STR="aprun -n %s " %(num_slave_nodes+1) 

		'''

		Notice: Your job was NOT submitted 

		Memory requests through PBS are not allowed on Kraken.
		A node's memory is accessible to all tasks running on the node.
		If a task requires all memory on a node, the PBS node/core request
		should be altered to request all cores on the required node(s). 

		Please remove the memory request (-l mem= or -l vmem=) and resubmit the job. 

	------------------------------------------------------------------------------------------------------------------

		The Kraken system is a Cray XT5 with 9,408 compute nodes interconnected with the SeaStar router through HyperTransport. The SeaStars are all interconnected in a 3-D torus topology. It is a massively parallel processing (MPP) machine. Each compute node has two six-core 2.6 GHz AMD Opterons for a total of 112,896 cores. All nodes have 16 Gbytes of DDR2 memory: 1.33 Gbytes of memory per core. There is 3.3 PB of storage. GRAM is deployed, but only as part of the science gateway kit. GRAM access must be preauthorized.

		Regular production work is placed in the small (0 to 512 cores), medium (513 to 8192 cores), and large (8193 to 49,536 cores) queues. Note, since queues are divided by size; one does not have to specify a queue type because it will be automatically determined. These queues have a 24 hour wall time.

		#!/bin/bash
		#PBS -A UT-NTNL0121
		#PBS -l size=192,walltime=01:35:00
		cd $PBS_O_WORKDIR
		aprun -n 192 ./a.out

		'''


	else:
		error_exit_with_message("Unsupported COMPUTER_CLUSTER_NAME=%s" %(COMPUTER_CLUSTER_NAME))

	print "ENVIROMENT_VARIABLE: COMPUTER_CLUSTER_NAME=%s | CORES_PER_NODE=%s | MPI_QUEUE_NAME=%s " %(COMPUTER_CLUSTER_NAME, CORES_PER_NODE, MPI_QUEUE_NAME)


	if(memory_reserve > MAX_MEMORY_PER_CORE):
		print "------------------------------------------------------------------------------------------------"
		print "WARNING: memory_reserve (%smb) > MAX_MEMORY_PER_CORE (%smb) on %s " %(memory_reserve, MAX_MEMORY_PER_CORE, COMPUTER_CLUSTER_NAME)
		print "WARNING: setting memory_reserve to MAX_MEMORY_PER_CORE!"
		memory_reserve=MAX_MEMORY_PER_CORE
		print "------------------------------------------------------------------------------------------------"


	#Submit the MPI job. This includes both the master_job and all the slave_jobs.
	qsub_submit_file="QSUB_MPI_DAG_JOB.sh" 

	if(exists(qsub_submit_file)): submit_subprocess("rm %s" %(qsub_submit_file))

	QSUB_JOB = open( qsub_submit_file, 'w')

	QSUB_JOB.write( '#!/bin/bash\n\n' )

	QSUB_JOB.write( '#PBS  -N %s\n\n' %(MPI_job_name))
	QSUB_JOB.write( '#PBS  -o %s\n\n' %(MPI_outfile + "_QSUB" ))
	QSUB_JOB.write( '#PBS  -e %s\n\n' %(MPI_errfile + "_QSUB" ))
	QSUB_JOB.write( '#PBS -l walltime=%d:00:00\n\n' %(wall_time))

	QSUB_JOB.write( '#PBS %s \n\n' %(PBS_REQUEST_COMPUTE_NODES_LINE) )   # number of nodes and number of processors per node requested. Would it be more efficient (shorter queue) time if don't don't demand 'whote' nodes?

	QSUB_JOB.write( '#PBS -M sripakpa@stanford.edu\n\n') # Address for email notification
	QSUB_JOB.write( '#PBS -m bea\n\n')  # "bea" mean mail is sent at the beginning, end, and at abort time (if it happens) 

	if(memory_reserve!=0):

		QSUB_JOB.write( '#PBS -l mem=%dmb\n\n' %(memory_reserve*(num_slave_nodes+1)))
		QSUB_JOB.write( '#PBS -l pmem=%dmb\n\n' %(memory_reserve)) #memory per CPU process. 

	if(MPI_QUEUE_NAME!=""):
		QSUB_JOB.write( '#PBS -q %s\n\n' %(MPI_QUEUE_NAME))

	QSUB_JOB.write( 'cd $PBS_O_WORKDIR\n\n' )

	QSUB_JOB.write( 'echo $JOB_ID > MPI_JOB_ID.txt\n\n')

	#If using OpenMPI, mpirun -np <procs> -hostfile $PBS_NODEFILE <executable>
	#If using mvapich2, mpirun_rsh -np 256 -hostfile $PBS_NODEFILE <executable>

	QSUB_JOB.write( '%s %s >%s 2>%s\n\n' %(MPI_EXEC_WRAPPER_STR, MPI_job_script, MPI_outfile, MPI_errfile)) #I think ibrun is specific command for submitting MPI jobs on the Lonestar XSEDE cluster.

	QSUB_JOB.close()

	submit_subprocess_allow_retry('qsub -V %s' %(qsub_submit_file)) #-V --> inherit the submission environment

