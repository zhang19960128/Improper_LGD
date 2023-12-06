#!/bin/bash
#SBATCH -n 96
#SBATCH -p all
#SBATCH -t 01:00:00
export ABI_PSPDIR=/workspace/jiahaoz/PACKAGE_INSTALL/abinit-9.4.1
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/workspace/jiahaoz/PACKAGE_INSTALL/abinit_support/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:workspace/jiahaoz/PACKAGE_INSTALL/abinit_support/libxc/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:workspace/jiahaoz/PACKAGE_INSTALL/curl/lib
cp /workspace/jiahaoz/PACKAGE_INSTALL/abinit/bin/abinit .
srun -n 36 -c2 --cpu_bind=cores ./abinit PATH7 > out
