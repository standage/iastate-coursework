#!/bin/tcsh

#PBS -o /home/standage/scratch/hw07output.txt
#PBS -e /home/standage/scratch/hw07error.txt
#PBS -lmem=256Mb,nodes=16:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
mpirun -np 8  hw07
mpirun -np 16 hw07
mpirun -np 32 hw07
