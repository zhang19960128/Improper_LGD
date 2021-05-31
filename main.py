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
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[Dw,Dv]=obtainmode.obtainmode(inputfiles.natomExp,masslist,inputfiles.dfptinExp,inputfiles.dfptoutExp);
ExpandPosition=analyzePH.readposition(axis,inputfiles.dfptinExp,inputfiles.natomExp);
symop=analyzePH.readsymmetry(inputfiles.dfptoutExp);
length=len(symop)
localsymop=symop[rank:length:size];
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
print('-------------------------GAMMA------------------------------------');
[primitiveGa,wGa,vGa]=primitive.obtainprimitivemode(inputfiles.natomPGa,inputfiles.dfptinPGa,inputfiles.dfptoutPGa,inputfiles.modePGaname);
EvGamma=extend.modemap(ExpandPosition,masslist,axis,primitiveGa,vGa.real,[1,1,1])
if rank==0:
  IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvGamma,inputfiles.modeExpGaname);
modematch.modecheck(EvGamma);
print(modematch.groupmode(wGa))
matchcoeffGamatrix=np.zeros((len(localsymop),len(EvGamma)));
for i in range(len(localsymop)):
  matchcoeffGa=modematch.modematch(localsymop[i],localTlist[i],EvGamma,wGa,axis);
  matchcoeffGamatrix[i]=np.copy(matchcoeffGa);
  startstr="";
  for j in range(len(matchcoeffGa)):
    startstr=startstr+" {:4.0f}".format(matchcoeffGamatrix[i][j])
  print("SymOP {0:15s}".format(SYMM.operations[i]),startstr)
print('-------------------------- M ------------------------------------');
[primitiveM,wM,vM]=primitive.obtainprimitivemode(inputfiles.natomPM,inputfiles.dfptinPM,inputfiles.dfptoutPM,inputfiles.modePMname);
EvM=extend.modemap(ExpandPosition,masslist,axis,primitiveM,vM.real,[-1,-1,-1]);
if rank==0:
  IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvM,inputfiles.modeExpMname);
modematch.modecheck(EvM);
print(modematch.groupmode(wM))
matchcoeffMmatrix=np.zeros((len(localsymop),len(EvM)));
for i in range(len(localsymop)):
  matchcoeffM=modematch.modematch(localsymop[i],localTlist[i],EvM,wM,axis,debug=0);
  matchcoeffMmatrix[i]=np.copy(matchcoeffM);
  startstr="";
  for j in range(len(matchcoeffM)):
    startstr=startstr+" {:4.0f}".format(matchcoeffMmatrix[i][j])
  print("SymOP {0:15s}".format(SYMM.operations[i]),startstr);
print('-------------------------- X ------------------------------------');
[primitiveX,wX,vX]=primitive.obtainprimitivemode(inputfiles.natomPX,inputfiles.dfptinPX,inputfiles.dfptoutPX,inputfiles.modePXname);
EvX=extend.modemap(ExpandPosition,masslist,axis,primitiveX,vX.real,[-1,-1,-1]);
if rank==0:
  IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvX,inputfiles.modeExpXname);
modematch.modecheck(EvX);
print(modematch.groupmode(wX))
matchcoeffXmatrix=np.zeros((len(localsymop),len(EvX)));
for i in range(len(localsymop)):
  matchcoeffX=modematch.modematch(localsymop[i],localTlist[i],EvX,wX,axis,debug=0);
  matchcoeffXmatrix[i]=np.copy(matchcoeffX);
  startstr="";
  for j in range(len(matchcoeffX)):
    startstr=startstr+" {:4.0f}".format(matchcoeffXmatrix[i][j])
  print("SymOP {0:15s}".format(SYMM.operations[i]),startstr);
print('--------------------------- LGD ----------------------------------');
print('--------------------------- First Order with X--------------------------');
for i in range(2,len(EvM)):
  diff=np.linalg.norm(matchcoeffMmatrix[:,j]-matchcoeffGamatrix[:,6]);
  if diff < 1e-6:
    print(i,diff)
#print(matchcoeffMmatrix[:,7]-matchcoeffGamatrix[:,6])
#print(matchcoeffMmatrix[:,7])
print('-------------------------- Second Order --------------------------');
for i in range(2,len(EvM)):
  for j in range(2,len(EvM)):
    diff=np.linalg.norm(matchcoeffMmatrix[:,i]*matchcoeffMmatrix[:,j]-matchcoeffGamatrix[:,6]);
    if diff < 1e-6:
      print(i,j);
