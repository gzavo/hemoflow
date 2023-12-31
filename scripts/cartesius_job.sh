#!/bin/bash
#SBATCH -t 0:30:00 #Ten minutes, adjust to need
#SBATCH -n 384    #number of needed processors, automatically calculate number of needed nodes
#SBATCH -p short #partition allocation
# #SBATCH --constraint=haswell # Runs only on Haswell nodes
		#name     max   type description
		# normal    5-00:00:00    720       thin      the default partition
		# fat       5-00:00:00     16       fat       fat node partition
		# short       01:00:00    720       thin      for short jobs
		# gpu       5-00:00:00     48       gpu       gpu nodes for production runs
		# gpu_short   01:00:00     64       gpu       gpu nodes for test runs
		# staging   5-00:00:00   1 core     server    for accessing the archive and external systems
#module rm mpi fortran c
#module load mpi/openmpi/1.10.2
#module load gcc/5.2.0
#Place the actual command here
#Dont forget to place this file in the directory you want to run from
time mpirun -n 384 ./hemoFlow setup.xml #dont forget to adjust -n !!
