import analyzePH
import numpy as np
import math
import cmath
import sys
import analyzePH
def exactmatch(position1,Poriginal,axis):
  natom=np.shape(position1)[0];
  exact=1;
  transform=np.zeros(natom);
  tick=0;
  for i in range(natom):
    mindist=1e10;
    for j in range(natom):
      dist=analyzePH.distance(position1[j],Poriginal[i],axis);
      if dist < mindist:
        mindist=dist;
        tick=j;
    if mindist < 1e-8:
      exact=exact*1;
      transform[i]=tick;
    else:
      exact=exact*0;
  return([exact,transform]);
def matchmass(transformlist,masslist):
  length=len(masslist);
  exact=1;
  for i in range(length):
    if np.abs(masslist[int(transformlist[i])]-masslist[i])<1e-6:
      exact=exact*1;
    else:
      exact=exact*0;
  return(exact);
def findexactmatch(OP,masslist,position,axis):
  sp=np.shape(position);
  PafterOP=np.zeros(sp);
  PafterOPandT=np.zeros(sp);
  transformlist=[];
  for m in range(sp[0]):
    temp=np.matmul(OP,position[m,0:3].reshape((3,1)));
    PafterOP[m]=np.copy(temp.reshape(3));
  for i in range(sp[0]):
    for j in range(sp[0]):
      dist=analyzePH.distance(PafterOP[i],position[j],axis);
      if dist < 1e-8 and np.abs(masslist[i]-masslist[j])<1e-8:
        transformlist.append(j);
        PafterOPandT[i]=np.copy(position[j]);
  if(len(transformlist)==sp[0]):
    pass;
  else:
    print("Cannot FIND MATCH")
  solution=[PafterOPandT,transformlist,PafterOP];
  return(solution)
def check(symop,axis,masslist,position):
  length=len(symop);
  sp=np.shape(position);
  transformlist=[];
  for i in range(length):
    OP=symop[i];
    [PafterOPandT,Transform,PafterOP]=findexactmatch(OP,masslist,position,axis);
    transformlist.append(Transform);
  return transformlist;
def translation(axis,masslist,qvector,position):
  sp=np.shape(position);
  PafterOP=np.zeros(sp);
  PafterOPandT=np.zeros(sp);
  transformlist=[];
  for m in range(sp[0]):
    for n in range(3):
      PafterOP[m]=PafterOP[m]+axis[n]*qvector[n];
    PafterOP[m]=PafterOP[m]+position[m];
  for i in range(sp[0]):
    for j in range(sp[0]):
      dist=analyzePH.distance(PafterOP[i],position[j],axis);
      if dist < 1e-8 and np.abs(masslist[i]-masslist[j])<1e-8:
        transformlist.append(j);
        PafterOPandT[i]=np.copy(position[j]);
  if(len(transformlist)==sp[0]):
    pass;
  else:
    print("Cannot FIND MATCH")
  switch=np.zeros((len(transformlist),len(transformlist)));
  for i in range(sp[0]):
    switch[transformlist[i]][i]=1;
  return([switch,transformlist]);
def modeafterOP(switch,natom,OP):
  Tmatrix=np.zeros((3*natom,3*natom));
  for i in range(natom):
    for j in range(natom):
      for m in range(3):
        for n in range(3):
          Tmatrix[3*i+m][3*j+n]=OP[m][n]*switch[i][j];
  return(Tmatrix);
def modeoperatingmatrix(symop,axis,masslist,qvector,position):
  Tlist=check(symop,axis,masslist,position);
  totalOP=[];
  natom=len(masslist);
  for t in range(len(qvector)):
    for i in range(len(symop)):
      trans=translation(axis,masslist,qvector[t],position)[0];
      transOP=modeafterOP(trans,natom,np.identity(3));
      switch=np.zeros((natom,natom));
      for k in range(natom):
        switch[Tlist[i][k]][k]=1;
      originalOP=modeafterOP(switch,natom,symop[i]);
      totalOP.append(np.matmul(transOP,originalOP));
  return totalOP;
if __name__=="__main__":
  pass;
else:
  pass;
