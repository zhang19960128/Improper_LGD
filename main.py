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
import sciconst
import modematch
import extend
import SYMM
import IOcharacter
#Note that symmetry operation also move the atoms#
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[Dw,Dv]=obtainmode.obtainmode(inputfiles.natomExp,masslist,inputfiles.dfptinExp,inputfiles.dfptoutExp);
ExpandPosition=analyzePH.readposition(axis,inputfiles.dfptinExp,inputfiles.natomExp);
symop=analyzePH.readsymmetry(inputfiles.dfptoutExp);
length=len(symop)
localsymop=symop[0:length];
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
print('-------------------------GAMMA------------------------------------');
[primitiveGa,wGa,vGa]=primitive.obtainprimitivemode(inputfiles.natomPGa,inputfiles.dfptinPGa,inputfiles.dfptoutPGa,inputfiles.modePGaname);
EvGamma=extend.modemap(ExpandPosition,masslist,axis,primitiveGa,vGa.real,[1,1,1])
IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvGamma,inputfiles.modeExpGaname);
modematch.modecheck(EvGamma);
print(modematch.groupmode(wGa))
matchcoeffGamatrix=np.zeros((len(localsymop),len(EvGamma)));
for i in range(len(localsymop)):
  matchcoeffGa=modematch.modematch(localsymop[i],localTlist[i],EvGamma,wGa,axis);
  matchcoeffGamatrix[i]=np.copy(matchcoeffGa);
IOcharacter.printoscreencharater(SYMM.operations,matchcoeffGamatrix)
print('-------------------------- M ------------------------------------');
[primitiveM,wM,vM]=primitive.obtainprimitivemode(inputfiles.natomPM,inputfiles.dfptinPM,inputfiles.dfptoutPM,inputfiles.modePMname);
EvM=extend.modemap(ExpandPosition,masslist,axis,primitiveM,vM.real,[-1,-1,-1]);
IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvM,inputfiles.modeExpMname);
modematch.modecheck(EvM);
print(modematch.groupmode(wM))
matchcoeffMmatrix=np.zeros((len(localsymop),len(EvM)));
for i in range(len(localsymop)):
  matchcoeffM=modematch.modematch(localsymop[i],localTlist[i],EvM,wM,axis);
  matchcoeffMmatrix[i]=np.copy(matchcoeffM);
IOcharacter.printoscreencharater(SYMM.operations,matchcoeffMmatrix)
FEmodeindex=6;
print('--------------------------- LGD Gamma*Gamma ----------------------------------');
for i in range(0,len(vGa)):
  diff=np.linalg.norm(matchcoeffGamatrix[:,i]-matchcoeffGamatrix[:,FEmodeindex]);
  if diff < 1e-6:
    print(i);
print('--------------------------- LGD Gamma*M -------------------------------------');
for i in range(0,len(vM)):
  diff=np.linalg.norm(matchcoeffMmatrix[:,i]-matchcoeffGamatrix[:,FEmodeindex]);
  if diff < 1e-6:
    print(i);
print('-------------------------- M*M --------------------------');
for i in range(0,len(EvM)):
  for j in range(i,len(EvM)):
    diff=np.linalg.norm(matchcoeffMmatrix[:,i]*matchcoeffMmatrix[:,j]-matchcoeffGamatrix[:,FEmodeindex]);
    if diff < 1e-6:
      print(i,j);
