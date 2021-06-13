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
import xlsxwriter
import pandas as pd
import ireducible
#import symmcharacter
#Note that symmetry operation also move the atoms#
symop=analyzePH.readsymmetry(inputfiles.dfptoutPGa);
length=len(symop)
localsymop=symop[0:length];
print("Number of Symmetry Operations:",len(localsymop))
[ExpandPosition,wmerge,vmerge]=primitive.obtainprimitivemode(inputfiles.natomExp,inputfiles.dfptinExp,inputfiles.dfptoutExp,inputfiles.modeExptotalname);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
print(modematch.groupmode(wmerge))
modematch.modecheck(vmerge);
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
matchcoeffmatrix=np.zeros((len(localsymop),len(vmerge)));
for i in range(0,len(localsymop)):
  matchcoeffGa=modematch.modematchsecond(localsymop[i],localTlist[i],vmerge,wmerge,axis);
  matchcoeffmatrix[i]=np.copy(matchcoeffGa);
datapd={"MOD"+str(i):np.round(matchcoeffmatrix[:,i],4) for i in range(len(vmerge))};
df = pd.DataFrame(data=datapd);
df.to_csv('Character.csv', index=False);
#--------- Mode*Mode -----#
modeindex=12;
for i in range(len(vmerge)):
  for j in range(i,len(vmerge)):
    if np.linalg.norm(matchcoeffmatrix[:,i].reshape(len(localsymop))*matchcoeffmatrix[:,j].reshape(len(localsymop))-matchcoeffmatrix[:,modeindex].reshape(len(localsymop))) < 1e-3:
      print(i,j)
''' Find 1D mode '''
[opmatrixinrep,opmatrix,irep]=ireducible.readirreducible();
print(irep)
print(opmatrix)
