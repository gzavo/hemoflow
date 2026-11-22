#!/bin/bash
# GPU build with several GPUs per node: gives every MPI rank of a node its own GPU, so that no rank
# keeps a CUDA context on another rank's GPU. Ranks beyond the number of GPUs share one (modulo);
# HemoFlow lets them wait. With Slurm, `srun --gpus-per-task=1` does the same without this script.
#
#   mpirun -np 8 scripts/gpu_bind.sh ./build-gpu/hemoFlow case.xml
#
# The number of GPUs per node comes from nvidia-smi, or from HEMOFLOW_GPUS_PER_NODE.
local=${OMPI_COMM_WORLD_LOCAL_RANK:-${MPI_LOCALRANKID:-${SLURM_LOCALID:-0}}}
gpus=${HEMOFLOW_GPUS_PER_NODE:-$(nvidia-smi -L 2>/dev/null | grep -c '^GPU')}
if [ "${gpus:-0}" -gt 0 ]; then
    export CUDA_VISIBLE_DEVICES=$((local % gpus))
fi
exec "$@"
