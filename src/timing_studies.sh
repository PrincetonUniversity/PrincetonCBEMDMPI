#!/bin/bash

tsteps=1000
n=216

for i in 1 2 
do
    time mpirun -np $i ./timing $tsteps LJ_$n.xml LJ.energy LJ_${n}_$i.out > timing_output_${n}_$i.out
    #time mpirun -np $i lmp_cluster -var timesteps $tsteps -var datafile LJ_$n.lammps -var outfile LJ_${n}_$i.out > timing_output_${n}_$i.out
done
    
