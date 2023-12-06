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
import constants
import printireducible
import positdiff
import modeproject
import plotprojection
#import symmcharacter
#Note that symmetry operation also move the atoms#
symop=analyzePH.readsymmetry(inputfiles.dfptoutExp);
length=len(symop)
localsymop=symop[0:length];
[ExpandPosition,wmerge,vmerge]=primitive.obtainprimitivemode(inputfiles.natomExp,inputfiles.dfptinExp,inputfiles.dfptoutExp,inputfiles.modeExptotalname);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
print('Axis is:=',axis)
natom=12;
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
Tlist=check.check(symop,axis,masslist,ExpandPosition);
#for i in range(len(Tlist)):
#  print(Tlist[i])
print(modematch.groupmode(wmerge))
modematch.modecheck(vmerge);
qvector=[
[0,0,0],
[0.0,0.5,0.5],
[0.5,0.0,0.5],
[0.5,0.5,0.0]
]
#print(check.translation(axis,masslist,qvector[0],ExpandPosition)[1])
#print(check.translation(axis,masslist,qvector[1],ExpandPosition))
#print(check.translation(axis,masslist,qvector[2],ExpandPosition))
#print(check.translation(axis,masslist,qvector[3],ExpandPosition))
modematrix=check.modeoperatingmatrix(localsymop,axis,masslist,qvector,ExpandPosition);
print("Length of symmetry is:=",len(modematrix))
printireducible.irrep(modematrix,vmerge,wmerge,axis);
replist=[0,2,3,4,5,6];
for i in replist:
  for j in replist:
    if i > j:
      continue;
    printireducible.mergeirrep(len(modematrix),wmerge,[i,j]);
printireducible.mergeirrep(len(modematrix),wmerge,[0,2,3,5]);
plotprojection.projectpath(natom,axis,masslist,vmerge)
modepro=modeproject.modeproject("./ABIOUTPUT/dfpt0666","./Pprofile/scf19",natom,axis,masslist,vmerge);
dp=positdiff.dp(natom,axis,"./ABIOUTPUT/dfpt0666","./Pprofile/scf19");
dp1D=dp.reshape(3*natom);
mass3D=[];
for i in range(natom):
  for j in range(3):
    mass3D.append(masslist[i]);
print(np.sqrt(np.array(mass3D)));
modelist=[2,11,13,15,19,21,22,29];
couplingtotal=[];
for i in modelist:
  for j in modelist:
    for k in modelist:
      for l in modelist:
        if i < j and j < k and k < l:
          couplingtotal.append([i,j,k,l]);
print(couplingtotal)
for i in couplingtotal:
  print("{0:4d} {1:4d} {2:4d} {3:4d}".format(i[0],i[1],i[2],i[3]))
print(len(couplingtotal))
scalecoup=[1,1,1,1,1];
for coup in couplingtotal:
  dp=np.zeros(np.shape(vmerge[0]));
  sampling=4;
  for t in range(4):
    if coup[t]==29:
      scalecoup[t]=-1.0;
    elif coup[t]==11:
      scalecoup[t]=-4.0;
    elif coup[t]==15:
      scalecoup[t]=-4.0;
    elif coup[t]==19:
      scalecoup[t]=-4.0;
    elif coup[t]==22:
      scalecoup[t]=-4.0;
    else:
      scalecoup[t]=4.0;
  for i in range(sampling+1):
    for j in range(sampling+1):
      for m in range(sampling+1):
        for n in range(sampling+1):
          dp=(i/(sampling+0.0)*vmerge[coup[0]]*scalecoup[0])+(j/(sampling+0.0)*vmerge[coup[1]]*scalecoup[1])+(m/(sampling+0.0)*vmerge[coup[2]]*scalecoup[2])+(n/(sampling+0.0)*vmerge[coup[3]]*scalecoup[3]);
          dp=dp/np.sqrt(np.array(mass3D));
          dp2D=dp.reshape(natom,3);
          positdiff.writeinterpolate(natom,axis,dp2D+ExpandPosition,"./phase/cubic","./MODEFIT4/"+"MODE"+str(coup[0])+"MODE"+str(coup[1])+"MODE"+str(coup[2])+"MODE"+str(coup[3])+"INTER"+str(i)+"INTER"+str(j)+"INTER"+str(m)+"INTER"+str(n));
