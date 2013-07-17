#!/bin/tcsh

#PBS -o /home/standage/scratch/hw08output.txt
#PBS -e /home/standage/scratch/hw08error.txt
#PBS -lmem=256Mb,nodes=16:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
mpirun -np 4  hw08
mpirun -np 8  hw08
mpirun -np 16 hw08
mpirun -np 32 hw08
