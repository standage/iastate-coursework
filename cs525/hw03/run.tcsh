#!/bin/tcsh

#PBS -o /home/standage/scratch/hw03output.txt
#PBS -e /home/standage/scratch/hw03error.txt
#PBS -lmem=512Mb,nodes=4:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
mpirun -np 8  hw03
