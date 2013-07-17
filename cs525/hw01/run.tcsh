#!/bin/tcsh

#PBS -o /home/standage/scratch/hw01output.txt
#PBS -e /home/standage/scratch/hw01error.txt
#PBS -lmem=256Mb,nodes=8:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
mpirun -np 4  hw01
mpirun -np 8  hw01
mpirun -np 16 hw01
