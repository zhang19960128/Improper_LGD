import numpy as np
import analyzePH
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
  w,v=np.linalg.eig(dyn2D);
  v=v.transpose();
  return [w,v];
