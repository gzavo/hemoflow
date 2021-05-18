#!/bin/bash
#SBATCH -t 0:30:00 # Time minutes, adjust to need
#SBATCH -n 384     # Number of needed processors, automatically calculate number of needed nodes
#SBATCH -p normal  # Partition allocation
#SBATCH --constraint=haswell # Runs only on Haswell nodes
		#name     max   type description
		# normal    5-00:00:00    720       thin      the default partition
		# fat       5-00:00:00     16       fat       fat node partition
		# short       01:00:00    720       thin      for short jobs
		# gpu       5-00:00:00     48       gpu       gpu nodes for production runs
		# gpu_short   01:00:00     64       gpu       gpu nodes for test runs
		# staging   5-00:00:00   1 core     server    for accessing the archive and external systems

# Sample job script to run on Cartesius (SURF).

BASE_DIR=$(pwd -P)
BASE_DATA="${BASE_DIR}/data"
BASE_EXT="*.xml"
BASE_EXE="${BASE_DIR}/hemoFlow"

echo $BASE_EXE 

#find $BASE_DATA -type f -name $BASE_EXT

for LOCAL_PATH in $(find $BASE_DATA -type f -name $BASE_EXT); do
    #echo $LOCAL_PATH
    LOCAL_DIR=${LOCAL_PATH%/*}
    LOCAL_FILE="${LOCAL_PATH##*/}"

    cd $LOCAL_DIR

	# Uncomment these to check the output before submitting the job
    #echo $LOCAL_DIR
    #echo $LOCAL_FILE

	# Execute the tasks
    mpirun -np 384 $BASE_EXE $LOCAL_FILE
    
done

#Move back to original root dir
cd $BASE_DIR
