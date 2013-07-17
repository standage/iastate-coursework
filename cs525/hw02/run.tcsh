#!/bin/tcsh

#PBS -o /home/standage/scratch/hw02output.txt
#PBS -e /home/standage/scratch/hw02error.txt
#PBS -lmem=256Mb,nodes=16:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
echo "!--------------------"
echo "! n=3500"
echo "!--------------------"
mpirun -np 4  hw02
mpirun -np 4  hw02
mpirun -np 4  hw02
mpirun -np 8  hw02
mpirun -np 8  hw02
mpirun -np 8  hw02
mpirun -np 16 hw02
mpirun -np 16 hw02
mpirun -np 16 hw02
mpirun -np 32 hw02
mpirun -np 32 hw02
mpirun -np 32 hw02
