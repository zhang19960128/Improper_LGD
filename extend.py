import numpy as np
import analyzePH
def gammamap(Hposition,masslist,axis):
  natom=len(Hposition);
  newaxis=np.copy(axis);
  for i in range(2):
    for j in range(2):
      newaxis[i][j]=axis[i][i]/2.0;
  newaxis[1][0]=-newaxis[1][0];
  pair=[];
  for i in range(natom):
    for j in range(i+1,natom):
      temp=Hposition[i]-Hposition[j];
      solution=np.matmul(np.linalg.inv(newaxis.transpose()),temp.reshape(3,1))
      match=1;
      for k in range(3):
        if np.abs(solution[k][0]-round(solution[k][0])) < 1e-6 :
          match=match*1;
        else:
          match=match*0;
      if match == 1 and np.abs(masslist[i]-masslist[j]) < 1e-6 :
        pair.append([i,j]);
  return pair
def modemap(ExpP,masslist,axis,PrimP,Pv,factor):
  pairs=gammamap(ExpP,masslist,axis);
  mapdict=dict();
  for i in range(len(pairs)):
    mapdict[pairs[i][0]]=pairs[i][1];
    mapdict[pairs[i][1]]=pairs[i][0];
  numofv=len(Pv);
  Ev=np.zeros((numofv,len(Pv[0])*2));
  maptoExp=[];
  for i in range(len(PrimP)):
    for j in range(len(ExpP)):
      dist=analyzePH.distance(PrimP[i],ExpP[j],axis);
      if dist < 1e-3:
        maptoExp.append(j);
      else:
        pass;
  for i in range(len(Ev)):
    for j in range(len(PrimP)):
      for k in range(3):
        Ev[i][maptoExp[j]*3+k]=Pv[i][3*j+k];
        Ev[i][mapdict[maptoExp[j]]*3+k]=Pv[i][3*j+k]*factor;
  for i in range(len(Ev)):
    Ev[i]=Ev[i]/np.linalg.norm(Ev[i]);
  return Ev;
