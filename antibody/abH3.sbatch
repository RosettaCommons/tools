#!/bin/bash
#----------------------------------------------------
# SLURM job script to run RosettaAntibody MPI application 
# on TACC's Stampede system.
#
# usage: from within the Ab directory, sbatch ../abH3.sbatch 
# (so cwd is the Ab directory)

# J. Gray - 2013 March 1
#----------------------------------------------------

#SBATCH -J ABNAME0000             # Job name
#SBATCH -o outerr/H3.%j.out       # Name of stdout output file (%j expands to jobId)
#SBATCH -e outerr/H3.%j.err       # Name of stdout output file (%j expands to jobId)
#SBATCH --mail-user=jeffreyjgray@gmail.com
#SBATCH --mail-type=FAIL
#SBATCH -A 454HTPSeq          # <-- Allocation name to charge job against

#development
###SBATCH -n 16                 # Total number of mpi tasks requested (16 cores/node)
###SBATCH -t 01:00:00           # Run time (hh:mm:ss)
###SBATCH -p development        # Queue name

#production
#SBATCH -n 64                 # Total number of mpi tasks requested (16 cores/node)
#SBATCH -t 24:00:00           # Run time (hh:mm:ss)
#SBATCH -p normal             # Queue name

ROSETTA=$WORK/git/Rosetta
ROSETTABIN=$ROSETTA/main/source/bin
ROSETTAEXE=antibody_H3
COMPILER=mklmpi.linuxiccrelease
EXE=$ROSETTABIN/$ROSETTAEXE.$COMPILER
echo Starting MPI job running $EXE

export MKL_MIC_ENABLE=1

time ibrun $EXE @../abH3.flags
