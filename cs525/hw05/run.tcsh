#!/bin/tcsh

#PBS -o /home/standage/scratch/hw05output.txt
#PBS -e /home/standage/scratch/hw05error.txt
#PBS -lmem=512Mb,nodes=16:ppn=2,cput=0:05:00,walltime=0:10:00

cd $PBS_O_WORKDIR
echo ""
echo "------ np=4 -----"
mpirun -np 4 hw05
echo ""
echo "------ np=8 -----"
mpirun -np 8 hw05
echo ""
echo "------ np=16 -----"
mpirun -np 16 hw05
echo ""
echo "------ np=32 -----"
mpirun -np 32 hw05
