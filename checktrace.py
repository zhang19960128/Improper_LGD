import numpy as np
from itertools import permutations
def matrixsame(IREPdat1,IREPdat2,numop):
  s=0.0;
  for i in range(numop):
    s=s+np.linalg.norm(IREPdat1[i]-IREPdat2[i]);
  return(s);
IREPdat=[];
symmopnum=192;
for i in range(9):
  IREP=np.loadtxt("IREP"+str(i)+".dat");
  shape=np.shape(IREP);
  xdim=int(shape[0]/shape[1]);
  symmopnum=xdim;
  ydim=int(shape[1]);
  IREPnew=IREP.reshape((xdim,ydim,ydim));
  IREPdat.append(IREPnew);
diff=np.zeros(symmopnum);
repone=0;
reptwo=4;
for i in range(symmopnum):
  diff[i]=np.trace(IREPdat[repone][i])-np.trace(IREPdat[reptwo][i]);
print(np.linalg.norm(diff))
'''
dim=6;
xdim=192;
ydim=6;
ydim=6;
permute= list(permutations(range(0,dim)));
changesignlist=[];
for i in range(-1,2,2):
  for j in range(-1,2,2):
    for k in range(-1,2,2):
      for m in range(-1,2,2):
        for n in range(-1,2,2):
          for l in range(-1,2,2):
            changesign=np.array([i,j,k,m,n,l]);
            changesignlist.append(changesign);
for i in range(len(permute)):
  permutation=permute[i];
  Tmatrix=np.zeros((dim,dim));
  MatAfterT=np.zeros((xdim,ydim,ydim));
  for m in range(len(changesignlist)):
    for j in range(dim):
      Tmatrix[permutation[j]][j]=1*changesignlist[m][j];
    for k in range(xdim):
      tempmatrix=np.matmul(Tmatrix.transpose(),IREPdat[reptwo][k]);
      tempmatrix=np.matmul(tempmatrix,Tmatrix);
      MatAfterT[k]=np.copy(tempmatrix);
    temp=matrixsame(IREPdat[repone],MatAfterT,xdim);
    if temp < 100:
      print(temp,Tmatrix)
'''
