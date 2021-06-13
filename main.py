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
import analyzeEXP
#Note that symmetry operation also move the atoms#
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[Dw,Dv]=obtainmode.obtainmode(inputfiles.natomExp,masslist,inputfiles.dfptinExp,inputfiles.dfptoutExp);
ExpandPosition=analyzePH.readposition(axis,inputfiles.dfptinExp,inputfiles.natomExp);
symop=analyzePH.readsymmetry(inputfiles.dfptoutPGa);
length=len(symop)
localsymop=symop[0:length];
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
print("Number of Symmetry Operations:",len(localsymop))
print('-------------------------GAMMA------------------------------------');
[primitiveGa,wGa,vGa]=primitive.obtainprimitivemode(inputfiles.natomPGa,inputfiles.dfptinPGa,inputfiles.dfptoutPGa,inputfiles.modePGaname);
EvGamma=extend.modemap(ExpandPosition,masslist,axis,primitiveGa,vGa.real,[1,1,1])
IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvGamma,inputfiles.modeExpGaname);
modematch.modecheck(EvGamma);
print(modematch.groupmode(wGa))
for i in range(len(wGa)):
  print("{0:10.7f}".format(wGa[i]))
matchcoeffGamatrix=np.zeros((len(localsymop),len(EvGamma)));
for i in range(0,len(localsymop)):
  matchcoeffGa=modematch.modematchsecond(localsymop[i],localTlist[i],EvGamma,wGa,axis);
  matchcoeffGamatrix[i]=np.copy(matchcoeffGa);
IOcharacter.printoscreencharater(SYMM.operations,matchcoeffGamatrix)
print('-------------------------- M ------------------------------------');
[primitiveM,wM,vM]=primitive.obtainprimitivemode(inputfiles.natomPM,inputfiles.dfptinPM,inputfiles.dfptoutPM,inputfiles.modePMname);
EvM=extend.modemap(ExpandPosition,masslist,axis,primitiveM,vM.real,[-1,-1,-1]);
IOmode.printmode(ExpandPosition,namelist,masslist,axis,EvM,inputfiles.modeExpMname);
modematch.modecheck(EvM);
print(modematch.groupmode(wM))
matchcoeffMmatrix=np.zeros((len(localsymop),len(EvM)));
for i in range(len(wM)):
  print("{0:10.7f}".format(wM[i]))
for i in range(len(localsymop)):
  matchcoeffM=modematch.modematchsecond(localsymop[i],localTlist[i],EvM,wM,axis);
  matchcoeffMmatrix[i]=np.copy(matchcoeffM);
IOcharacter.printoscreencharater(SYMM.operations,matchcoeffMmatrix)
FEmodeindex=8;
for i in range(12):
  print("{0:10.7f} {1:10.7f}".format(wGa[i],wM[i]))
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
McoupleM=[];
for i in range(0,len(EvM)):
  for j in range(i,len(EvM)):
    diff=np.linalg.norm(matchcoeffMmatrix[:,i]*matchcoeffMmatrix[:,j]-matchcoeffGamatrix[:,FEmodeindex]);
    if diff < 1e-6:
      print(i,j);
      McoupleM.append([i,j]);
'''
[pTCtoGa,pTCtoM,pOCtoGa,pOCtoM]=analyzeEXP.analyzetransition(EvGamma,EvM);
IOmode.printvector(pTCtoGa);
IOmode.printvector(pTCtoM);
IOmode.printvector(pOCtoGa);
IOmode.printvector(pOCtoM);
'''
