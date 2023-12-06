import numpy as np
import IOmode
import analyzePH
import obtainmode
import sequence
import cmath
import sciconst
import modematch
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
def reorthogonalize(wmode,vmode):
  gplist=modematch.groupmode(wmode);
  for i in range(len(gplist)):
    if i==1:
      continue;
    print("Group ID is:=",i)
    if len(gplist[i])==3:
      natom=len(vmode[0])/3;
      natom=int(natom);
      newmode=[];
      for j in range(3): # j is the direction we want to normalize;
        total=[0,1,2];
        total.remove(j);
        v1=np.zeros(2*natom);
        v2=np.zeros(2*natom);
        v3=np.zeros(2*natom);
        v1[0:natom]=np.copy(vmode[gplist[i][0]][total[0]:3*natom:3]);
        v1[(natom):(2*natom)]=np.copy(vmode[gplist[i][0]][total[1]:3*natom:3]);
        v2[0:natom]=np.copy(vmode[gplist[i][1]][total[0]:3*natom:3]);
        v2[(natom):(2*natom)]=np.copy(vmode[gplist[i][1]][total[1]:3*natom:3]);
        v3[0:natom]=np.copy(vmode[gplist[i][2]][total[0]:3*natom:3]);
        v3[(natom):(2*natom)]=np.copy(vmode[gplist[i][2]][total[1]:3*natom:3]);
        #solve C1*V1+C2*V2+C3*V3=0, let C1=1.0;
        Dmatrix=np.zeros((2,2*natom));
        Dmatrix[0]=np.copy(v2);
        Dmatrix[1]=np.copy(v3);
        C1=1.0
        A=np.matmul(Dmatrix,Dmatrix.transpose());
        b=-C1*np.matmul(np.copy(v1),Dmatrix.transpose());
        result=np.matmul(b,np.linalg.inv(A));
        vnew=vmode[gplist[i][0]]*C1+vmode[gplist[i][1]]*result[0]+vmode[gplist[i][2]]*result[1];
        newmode.append(vnew);
      for j in range(3):
        vmode[gplist[i][j]]=np.copy(newmode[j]/np.linalg.norm(newmode[j]));
def obtainprimitivemode(natom,dfptin,dfptout,modename):
  [masslist,namelist]=sequence.sequence(natom);
  [w,v]=obtainmode.obtainmode(natom,masslist,dfptin,dfptout);
  reorthogonalize(w,v);
  path='./TOTAL/basishalf';
  for i in range(18,24):
    data=np.loadtxt(path+'/MODE'+str(i)+".back");
    sp=np.shape(data);
    newflat=data.reshape(sp[0]*sp[1]);
    vtemp=np.copy(newflat/np.linalg.norm(newflat));
    sumtemp=0;
    for j in range(18,24):
      sumtemp=sumtemp+(vtemp.dot(v[j]))**2;
    if np.abs(sumtemp-1.0) < 1e-5:
      print("Manual Assign is good!!!");
  for i in range(18,24):
    data=np.loadtxt(path+'/MODE'+str(i)+'.back');
    sp=np.shape(data);
    newflat=data.reshape(sp[0]*sp[1]);
    v[i]=np.copy(newflat/np.linalg.norm(newflat));
  for i in range(6,12):
    data=np.loadtxt(path+'/MODE'+str(i)+'.back');
    sp=np.shape(data);
    newflat=data.reshape(sp[0]*sp[1]);
    vtemp=np.copy(newflat/np.linalg.norm(newflat));
    sumtemp=0;
    for j in range(6,12):
      sumtemp=sumtemp+(vtemp.dot(v[j]))**2;
    if np.abs(sumtemp-1.0) < 1e-5:
      print("Manual assign is good!!!")
  for i in range(6,12):
    data=np.loadtxt(path+'/MODE'+str(i)+'.back');
    sp=np.shape(data);
    newflat=data.reshape(sp[0]*sp[1]);
    v[i]=np.copy(newflat/np.linalg.norm(newflat));
  for i in range(24,30):
    data=np.loadtxt(path+'/MODE'+str(i)+'.back');
    sp=np.shape(data);
    newflat=data.reshape(sp[0]*sp[1]);
    vtemp=np.copy(newflat/np.linalg.norm(newflat));
    sumtemp=0;
    for j in range(24,30):
      sumtemp=sumtemp+(vtemp.dot(v[j]))**2;
    if np.abs(sumtemp-1.0) < 1e-5:
      print("Manual assign is good!!!")
  for i in range(24,30):
    data=np.loadtxt(path+'/MODE'+str(i)+'.back');
    sp=np.shape(data);
    newflat=data.reshape(sp[0]*sp[1]);
    v[i]=np.copy(newflat/np.linalg.norm(newflat));
  axis=analyzePH.readaxis(dfptin);
  Hposition=analyzePH.readposition(axis,dfptin,natom);
  IOmode.printmode(Hposition,namelist,masslist,axis,v,modename);
#  rotationmatrix=np.array([[np.sqrt(2)/2,np.sqrt(2)/2,0],[-np.sqrt(2)/2,np.sqrt(2)/2,0],[0,0,1]]);
#  scale=np.array([np.sqrt(2),np.sqrt(2),1]);
#  [primitive,axis]=rotate(Hposition,axis,rotationmatrix,scale);
  primitive=np.copy(Hposition);
  return [primitive,w,v];
