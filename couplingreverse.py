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
modelist=[11,13,15,19,21,22,29];
modeproject=modepro[0].copy();
dp=np.zeros(3*natom);
mass3D=[];
for i in range(natom):
  for j in range(3):
    mass3D.append(masslist[i]);
print(np.sqrt(np.array(mass3D)));
for i in range(3*natom):
  dp=dp+modeproject[i]*vmerge[i];
dp=dp/np.sqrt(np.array(mass3D));
dp2D=dp.reshape(natom,3);
positdiff.writeinterpolate(natom,axis,dp2D+ExpandPosition,"./phase/cubic","reverse.abinit");
