#!/bin/bash
#SBATCH -t 15-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1 
#SBATCH --cpus-per-task=40
#SBATCH --threads-per-core=2
#SBATCH -e slurm-%j.err
#SBATCH -o slurm-%j.out
##SBATCH --exclusive
##SBATCH --nodelist=cl9

# this example will allocate 40 processors in a node to be fully used by one multithreading matlab instance

# run with:
# sbatch slurm_submission.sh

echo $HOME

matlab -nosplash -nodisplay -batch "pid=0;call_main" -logfile rift2ridge2D.log
#sleep 5
#srun -N 1 python $HOME/docs/programming/hpc/example_serial/demo.py 1
