#!/bin/bash
#SBATCH -n 96
#SBATCH -q debug
#SBATCH -t 00:60:00
export ABI_PSPDIR=/workspace/jiahaoz/PACKAGE_INSTALL/abinit-9.4.1 #/JTH-PBE-atomicdata-1.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/workspace/jiahaoz/PACKAGE_INSTALL/abinit_support/lib
abinitexec=/workspace/jiahaoz/PACKAGE_INSTALL/abinit/bin/
cp $abinitexec/abinit 
/workspace/jiahaoz/PACKAGE_INSTALL/INTEL/mpi/2021.2.0/bin/mpirun -n 28 ./abinit MODE19MODE21MODE22INTER6INTER6INTER6
