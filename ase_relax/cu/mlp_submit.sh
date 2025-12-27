#!/bin/bash
#SBATCH --account=rrg-ovoznyy
#SBATCH --job-name=md_ase
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=192
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
##SBATCH --mail-type=END,FAIL,ALL
##SBATCH --mail-user=hgabbas71@gmail.com
#SBATCH --output=job_output.log
#SBATCH --error=job_error.log

module load gcc/12.3
module load openmpi/4.1.5 
module load mpi4py/3.1.4
source /home/abbas71/ase/bin/activate
export PYTHONPATH=/home/abbas71/ase/bin:$PYTHONPATH
export UCX_TLS=tcp
# FIX: Disable the problematic UCX feature
export UCX_VFS_ENABLE=n
export UCX_MEM_HOOKS=no
python3 relax.py > output.log 2>&1
