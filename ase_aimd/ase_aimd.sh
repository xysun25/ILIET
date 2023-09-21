#!/bin/bash  
#SBATCH -p amd_256
#SBATCH -N 1              
#SBATCH -n 64
#SBATCH -J ase_aimd

# For vasp
module load intel/18.0.2-thc mpi/intel/18.0.2-thc-public4
export PATH="/public4/home/sc57227/xysun/software/vasp/5.4.4/bin/":$PATH
ulimit -s unlimited
NP=$[$SLURM_NPROCS * $SLURM_JOB_NUM_NODES]

# For ase
export ASE_VASP_COMMAND="mpirun -np $NP vasp_gam"
export VASP_PP_PATH="/public4/home/sc57227/xysun/software/vasp/5.4.4/potential"

python -u ase_aimd.py
