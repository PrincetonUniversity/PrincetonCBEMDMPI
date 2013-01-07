FPAPC5242012
============

Final Project for APC 524

README for CBEMD: Parallelized Molecular Dynamics in Various Thermodynamic Ensembles

Dependencies: Boost, Google Tests, OpenMPI

Boost---------------------------------------------------------------------------
Boost must be downloaded and installed locally in order for the code to operate. 
You can download Boost from http://www.boost.org/

Once installed, in Makefile and MakefileTests, you must point the variable 

PATHTOBOOST to where boost is located.

For example, if boost is located in /home/gkhoury/boost_1_52_0

set PATHTOBOOST = /home/gkhoury/boost_1_52_0

in both Makefile and MakefileTests


Google Tests--------------------------------------------------------------------
Google Tests can be downloaded from 
http://code.google.com/p/googletest/downloads/list
Our code was tested using version 1.6.0. It must be installed prior to making 
MakefileTests. Information about installation can be found
on the google test webpage.

Once installed, in MakefileTests, you must point the variable 

GTESTDIR to the path where googletests are located. 

For example, if googletests is located in /home/gkhoury/gtest-1.6.0

set GTESTDIR = /home/gkhoury/gtest-1.6.0

in MakefileTests


OpenMPI-------------------------------------------------------------------------
For successful compilation and running, we recommend compiling using

openmpi/gcc/1.2.8/64
which can be downloaded at http://www.open-mpi.org/software/ompi/v1.2/

and is already available on the university clusters sesame, adroit, and della.

Installation--------------------------------------------------------------------

cd to /src/

type make

Will compile ./verlet ./andersen 

type make -f MakefileTests all

Will compile ./tests

Execution-----------------------------------------------------------------------

Integraters
Verlet (NVE) Usage
./verlet nsteps dt xml_file energy_file output_file

Example Usage:
mpiexec -np 4 ./verlet 10 0.0005 LJ_1000.xml LJ.energy output

Andersen (NVT) Usage
./andersen nsteps dt xml_file energy_file output_file temperature nu

mpirun -n 4 ./andersen 10 0.0005 LJ_1000.xml LJ.energy outputandersen 1 10

to clean, type make clean

Note, the code was also tested by compiling with
openmpi/gcc/1.3.3/64
openmpi/intel-10.1/1.2.8/64
openmpi/intel-10.1/1.3.3/64
openmpi/intel-11.1/1.2.8/64
openmpi/intel-11.1/1.3.3/64
openmpi/intel-11.1/1.4.3/64

but is not gauranteed to work with them. To be safe, we recommend 
using openmpi/gcc/1.2.8/64


Sample PBS submission script for Della, Sesame, Tiger---------------------
-1. follow instructions above to link makefile to googletests and boost
0. type   module load openmpi 
1. cd to src and type make
2. type make -f MakefileTests all
3. copy the below into a file run.sh, edit parameters to liking
4. chmod u+x run.sh
5. qsub run.sh

#!/bin/bash -x
#PBS -r n
#PBS -l nodes=1:ppn=2,walltime=00:11:10

cd $PBS_O_WORKDIR

module load openmpi

mpiexec -np 2 ./andersen 10 0.0005 LJ_1000.xml LJ.energy outputandersen 10 50

mpiexec -np 2 ./verlet 10 0.0005 LJ_1000.xml LJ.energy output

