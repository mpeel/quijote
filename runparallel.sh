# #!/bin/bash
# #PBS -l nodes=1:ppn=2
# #PBS -l walltime=3:00:00:00
# #PBS -N cbass_spawn_reduction
# #PBS -m abe
# #PBS -M michael.peel@manchester.ac.uk

# source ~/.bashrc

# cd $PBS_O_WORKDIR

# echo "Running job in"
# echo `pwd`
# echo "List of nodes being used";
# cat $PBS_NODEFILE

# export CBASS_HOME="$PBS_O_WORKDIR/../"

# mpirun -np 2 --bynode --mca mpi_yield_when_idle 1 python run_q_smoothnoisemaps_parallel.py
mpirun -np 4 --mca mpi_yield_when_idle 1 python run_q_smoothnoisemaps_parallel.py
