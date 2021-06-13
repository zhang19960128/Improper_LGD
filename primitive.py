import numpy as np
import IOmode
import analyzePH
import obtainmode
import sequence
import cmath
import sciconst
def rotate(Hposition,axis,operatingmatrix,scale):
  newaxis=np.copy(axis);
  for i in range(3):
    newaxis[i]=np.matmul(operatingmatrix,axis[i].reshape((3,1))).reshape(3);
  crystalposition=np.copy(Hposition);
  newposition=np.copy(Hposition);
  natom=len(Hposition);
  for i in range(natom):
    crystalposition[i]=np.matmul(Hposition[i],np.linalg.inv(newaxis));
  for i in range(3):
    for j in range(3):
      if i!=j:
        newaxis[i][j]=0.0;
      else:
        newaxis[i][i]=np.linalg.norm(axis[i]);
  for i in range(natom):
    newposition[i]=np.matmul(crystalposition[i],newaxis);
  for i in range(3):
    newaxis[i][i]=newaxis[i][i]*scale[i];
  for i in range(natom):
    for j in range(3):
      if newposition[i][j]<0:
        newposition[i][j]=newposition[i][j]+newaxis[j][j];
      elif newposition[j][j]>newaxis[j][j]:
        newposition[i][j]=newposition[i][j]-newaxis[j][j];
  return [newposition,newaxis];
def obtainprimitivemode(natom,dfptin,dfptout,modename):
  [masslist,namelist]=sequence.sequence(natom);
  [w,v]=obtainmode.obtainmode(natom,masslist,dfptin,dfptout);
  axis=analyzePH.readaxis(dfptin);
  Hposition=analyzePH.readposition(axis,dfptin,natom);
  IOmode.printmode(Hposition,namelist,masslist,axis,v,modename);
#  rotationmatrix=np.array([[np.sqrt(2)/2,np.sqrt(2)/2,0],[-np.sqrt(2)/2,np.sqrt(2)/2,0],[0,0,1]]);
#  scale=np.array([np.sqrt(2),np.sqrt(2),1]);
#  [primitive,axis]=rotate(Hposition,axis,rotationmatrix,scale);
  primitive=np.copy(Hposition);
  return [primitive,w,v];
