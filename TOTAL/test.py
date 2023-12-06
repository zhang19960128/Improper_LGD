import numpy as np
modepaper=np.loadtxt("MODE23.back");
dm=np.shape(modepaper);
length=dm[0]*dm[1];
modepaper=modepaper.reshape(length);
modepaper=modepaper/np.linalg.norm(modepaper);
dim=3;
start=33;
mat1=np.zeros((dim,length));
#for i in range(dim):
#  vectorlist[i]=np.loadtxt("TOTAL"+str(int(start+i)));
#  mat1[i]=np.copy(vectorlist[i]);
A=np.matmul(mat1,mat1.transpose());
b=np.matmul(modepaper,mat1.transpose());
vectorlist=[];
for i in range(36):
  dat=np.loadtxt("TOTAL"+str(int(i)));
  vectorlist.append(dat);
for i in range(36):
  print(vectorlist[i].dot(modepaper))
#print(modepaper)
#print(A)
#print(np.matmul(b,np.linalg.inv(A)));
