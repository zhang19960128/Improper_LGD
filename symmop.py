import numpy as np
def symoperate(OP,modeindex,v,axis):
  sp=np.shape(v);
  VafterOP=np.zeros(np.shape(v));
  for i in range(int(sp[0]/3)):
    temp=np.matmul(OP,v[3*i:(3*i+3)].reshape((3,1)));
    for k in range(3):
      VafterOP[3*i+k]=temp[k,0];
  VafterP=np.zeros(np.shape(v));
  for i in range(int(sp[0]/3)):
    for k in range(3):
      VafterP[3*int(modeindex[i])+k]=VafterOP[3*i+k];
  return VafterP;
def multipleOperate(v,OPlist):
  length=len(OPlist);
  re=symoperate(OPlist[0],v);
  for i in range(1,length):
    re=symoperate(OPlist[i],re);
  return re;
def symoperate_matrix_op(OP,modeindex,v):
  natom=len(modeindex);
  switch=np.zeros((natom,natom));
  for i in range(len(modeindex)):
    switch[i][modeindex[i]]=1;
  Tmatrix=np.zeros((3*natom,3*natom));
  for i in range(natom):
    for j in range(natom):
      for m in range(3):
        for n in range(3):
          Tmatrix[3*i+m][3*j+n]=OP[m][n]*switch[i][j];
  return(np.matmul(Tmatrix,v.reshape(len(v),1)).reshape(len(v)));
