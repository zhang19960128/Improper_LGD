import numpy as np
import analyzePH
import modematch
def obtainmode(natom,masslist,dfptin,dfptout):
  dyn=analyzePH.obtaindyn(natom,dfptout);
  for i in range(natom):
    for j in range(natom):
      for m in range(3):
        for n in range(3):
          dyn[i][j][m][n]=dyn[i][j][m][n]/np.sqrt(masslist[i]*masslist[j]);
  dyn2D=np.zeros((3*natom,3*natom));
  for i in range(natom):
    for j in range(natom):
      dyn2D[(i*3):((i+1)*3),(j*3):((j+1)*3)]=np.copy(dyn[i,j,0:3,0:3]);
  w,v=np.linalg.eigh(dyn2D);
  v=v.transpose();
  v=schmidt(w,v);
  return [w,v];
def schmidt(w,v):
  wgroup=modematch.groupmode(w);#Gram-Schmidt Process#
  re=np.copy(v);
  for i in range(len(wgroup)):
    if len(wgroup[i])==1:
      pass;
    else:
      group=wgroup[i];
      for j in range(len(group)):
        for k in range(0,j):
          re[group[j]]=re[group[j]]-re[group[k]].dot(re[group[j]])*re[group[k]];
        re[group[j]]=np.copy(re[group[j]]/np.linalg.norm(re[group[j]]));
  return re;
