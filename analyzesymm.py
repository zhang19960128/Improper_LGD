import analyzePH
import numpy as np
import math
import cmath
import sys
import check
import IOmode
import inputfiles
from mpi4py import MPI
Ha=2*13.59*1.602176634*10**(-19);
bohr=0.529177*10**-10;
aumass=1.660539066*10**-27;
#Note that symmetry operation also move the atoms#
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
def modematch(v1,v2):
  length=len(v1);
  se=0.0;
  for i in range(length):
    se=se+v1[i]*v2[i];
  return se;
def findmodematch(v1,v):
  result=[];
  n=len(v);
  for i in range(n):
    re=modematch(v1,v[i]);
    if re**2 >0.8:
      result.append(i);
  return(result);
def groupmode(w):
  length=len(w);
  modeset=[i for i in range(length)];
  groupset=[];
  groupset.append([0]);
  modeset.remove(0);
  while(len(modeset)>0):
    length=len(groupset);
    tick=-1;
    for i in range(length):
      if np.abs(w[groupset[i][0]]-w[modeset[0]])/np.abs(w[groupset[i][0]])< 5*1e-2:
        groupset[i].append(modeset[0]);
        tick=tick+1;
    if np.abs(tick+1)<1e-6:
      groupset.append([modeset[0]]);
    modeset.remove(modeset[0]);
  return(groupset)
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
def sequence(natom):
  masslist=[];
  namelist=[];
  for i in range(int(1/6*natom)):
    masslist.append(178.49);
    namelist.append("Hf")
  for i in range(int(1/6*natom),int(2/6*natom)):
    masslist.append(91.224);
    namelist.append("Zr")
  for i in range(int(2/6*natom),natom):
    masslist.append(15.9999);
    namelist.append("O")
  return([masslist,namelist]);
def obtainprimitivemode(natom,dfptin,dfptout,modename):
  [masslist,namelist]=sequence(natom);
  [w,v]=obtainmode(natom,masslist,dfptin,dfptout);
  axis=analyzePH.readaxis(dfptin);
  Hposition=analyzePH.readposition(axis,dfptin,natom);
  if rank==0:
    for i in range(len(w)):
      IOmode.printmode(Hposition,namelist,axis,[v[i]],modename+"{0:6.2f}".format(cmath.sqrt(w[i])*np.sqrt(Ha/bohr/bohr/aumass)*10**-12*33.35641/2/3.141592653));
  rotationmatrix=np.array([[np.sqrt(2)/2,np.sqrt(2)/2,0],[-np.sqrt(2)/2,np.sqrt(2)/2,0],[0,0,1]]);
  scale=np.array([np.sqrt(2),np.sqrt(2),1]);
  [primitive,axis]=rotate(Hposition,axis,rotationmatrix,scale);
  return [primitive,w,v];
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
def findmodematchsubspace(v1,v,subgroup):
  mat=np.zeros((len(subgroup),len(v1)));
  for i in range(len(subgroup)):
    mat[i]=np.copy(v[subgroup[i]]);
  A=np.matmul(mat,mat.transpose());
  b=np.matmul(v1,mat.transpose());
  coeff=np.matmul(b,np.linalg.inv(A));
  subv=np.zeros(v1.shape);
  for i in range(len(coeff)):
    subv=subv+coeff[i]*v[subgroup[i]];
  subv=subv/np.linalg.norm(subv);
  return([subv.dot(v1),subv])
def groupmatchoperation(w,symop,Tlist,Ev,axis):
  gp=groupmode(w);
  matchcoeff=[];
  matchmode=[];
  for i in range(len(gp)):
    groupindex=i;
    vOP=symoperate(symop,Tlist,Ev[gp[groupindex][0]],axis);
    re=findmodematchsubspace(vOP,Ev,gp[groupindex]);
    for j in range(20):
      vOPafter=symoperate(symop,Tlist,re[1],axis);
      re=findmodematchsubspace(vOPafter,Ev,gp[groupindex]);
    if np.abs(re[0]) > 0.80:
      matchcoeff.append(re[0]);
      matchmode.append(re[1]);
    else:
      pass;
  return([matchcoeff,matchmode])
comm=MPI.COMM_WORLD
size=comm.Get_size();
rank=comm.Get_rank();
[primitiveM,wM,vM]=obtainprimitivemode(inputfiles.natomPM,inputfiles.dfptinPM,inputfiles.dfptoutPM,inputfiles.modePMname);
[primitiveGa,wGa,vGa]=obtainprimitivemode(inputfiles.natomPGa,inputfiles.dfptinGa,inputfiles.dfptoutGa,inputfiles.modePGaname);
[masslist,namelist]=sequence(inputfiles.natomExp);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[Dw,Dv]=obtainmode(inputfiles.natomExp,masslist,inputfiles.dfptinExp,inputfiles.dfptoutExp);
ExpandPosition=analyzePH.readposition(axis,inputfiles.dfptinExp,inputfiles.natomExp);
EvM=modemap(ExpandPosition,masslist,axis,primitiveM,vM.real,-1);
EvGamma=modemap(ExpandPosition,masslist,axis,primitiveGa,vGa.real,1)
if rank==0:
  IOmode.printmode(ExpandPosition,namelist,axis,EvM,inputfiles.modeExpMname);
  IOmode.printmode(ExpandPosition,namelist,axis,EvGamma,inputfiles.modeExpGaname);
  for i in range(len(wGa)):
    IOmode.printmode(ExpandPosition,namelist,axis,[EvGamma[i]],inputfiles.modeExpGaname+"{0:6.2f}".format(cmath.sqrt(wGa[i])*np.sqrt(Ha/bohr/bohr/aumass)*10**-12*33.35641/2/3.141592653));
symop=analyzePH.readsymmetry(inputfiles.dfptoutExp);
length=len(symop)
localsymop=symop[rank:length:size];
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
print(localsymop[3])
[matchcoeff,matchmode]=groupmatchoperation(wGa,localsymop[3],localTlist[3],EvGamma,axis);
IOmode.printmode(ExpandPosition,namelist,axis,matchmode,"./MATCH/MATCHMODE");
print(matchcoeff)
