import analyzePH
import numpy as np
import math
import cmath
import sys
import check
import IOmode
import inputfiles
import primitive
import sequence
import obtainmode
from mpi4py import MPI
import sciconst
import modematch
import extend
import SYMM
#Note that symmetry operation also move the atoms#
comm=MPI.COMM_WORLD
size=comm.Get_size();
rank=comm.Get_rank();
[primitiveM,wM,vM]=primitive.obtainprimitivemode(inputfiles.natomPM,inputfiles.dfptinPM,inputfiles.dfptoutPM,inputfiles.modePMname);
[primitiveGa,wGa,vGa]=primitive.obtainprimitivemode(inputfiles.natomPGa,inputfiles.dfptinGa,inputfiles.dfptoutGa,inputfiles.modePGaname);
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[Dw,Dv]=obtainmode.obtainmode(inputfiles.natomExp,masslist,inputfiles.dfptinExp,inputfiles.dfptoutExp);
ExpandPosition=analyzePH.readposition(axis,inputfiles.dfptinExp,inputfiles.natomExp);
EvM=extend.modemap(ExpandPosition,masslist,axis,primitiveM,vM.real,-1);
EvGamma=extend.modemap(ExpandPosition,masslist,axis,primitiveGa,vGa.real,1)
if rank==0:
  IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvM,inputfiles.modeExpMname);
  IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvGamma,inputfiles.modeExpGaname);
symop=analyzePH.readsymmetry(inputfiles.dfptoutExp);
length=len(symop)
localsymop=symop[rank:length:size];
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
modematch.modecheck(EvGamma);
modematch.modecheck(EvM);
print(modematch.groupmode(wGa))
print('-------------------------GAMMA------------------------------------');
for i in range(len(localsymop)):
  matchcoeff=modematch.modematch(localsymop[i],localTlist[i],EvGamma,wGa,axis);
  startstr="";
  for j in range(len(matchcoeff)):
    startstr=startstr+" {:4.0f}".format(matchcoeff[j])
  print("SymOP {0:15s}".format(SYMM.operations[i]),startstr)
print('-------------------------  M  ------------------------------------');
print(modematch.groupmode(wM))
for i in range(len(localsymop)):
  matchcoeff=modematch.modematch(localsymop[i],localTlist[i],EvM,wM,axis,debug=0);
  startstr="";
  for j in range(len(matchcoeff)):
    startstr=startstr+" {:4.0f}".format(matchcoeff[j])
  print("SymOP {0:15s}".format(SYMM.operations[i]),startstr)
