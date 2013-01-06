#! /bin/bash
#PBS -l nodes=1:ppn=1,walltime=1:00:00
#PBS -q caf

cd $PBS_O_WORKDIR


#source /home/gkhoury/.bashrc

module purge

module load openmpi/gcc/1.2.8/64

mpirun -n 1 ./andersen 10 0.0005 LJ_1000.xml LJ.energy outputandersen 10

#mpiexec -np 1 ./cbemd 10 0.0005 LJ_1000.xml LJ.energy output


