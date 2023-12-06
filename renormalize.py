import numpy as np
Tmatrix=np.loadtxt("Tmatrix.dat");
path='/workspace/jiahaoz/HZO/LG_Model/GroundP/HfO4/WaveVector/process/Improper_LGD/TOTAL/basishalf';
basis=[];
for i in range(24,30):
  data=np.loadtxt(path+'/MODE'+str(i)+'.back');
  sp=np.shape(data);
  newflat=data.reshape(sp[0]*sp[1]);
  vtemp=np.copy(newflat/np.linalg.norm(newflat));
  basis.append(vtemp);
newbasis=[];
for i in range(6):
  bas=np.zeros(6);
  bas[i]=1.0;
  bastwo=bas.reshape((6,1));
  result=np.matmul(Tmatrix,bastwo);
  basnow=result.reshape(6);
  result=np.zeros(36);
  for k in range(6):
    result=result+basnow[k]*basis[k];
  f=open(path+'/newbasis'+str(i)+'.dat','w');
  for k in range(12):
    f.write(str(result[3*k+0])+" "+str(result[3*k+1])+" "+str(result[3*k+2])+"\n");
  f.close();
