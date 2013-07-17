#!/bin/tcsh

#PBS -o /home/standage/scratch/hw06output.txt
#PBS -e /home/standage/scratch/hw06error.txt
#PBS -lmem=256Mb,nodes=16:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
echo "!--------------------"
echo "! n=3500"
echo "!--------------------"
mpirun -np 4  hw06
mpirun -np 8  hw06
mpirun -np 16 hw06
mpirun -np 32 hw06
