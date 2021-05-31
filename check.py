import analyzePH
import numpy as np
import math
import cmath
import sys
import analyzePH
def check(symop,axis,position):
  length=len(symop);
  sp=np.shape(position);
  for i in range(length):
    OP=symop[i];
    PafterOP=np.zeros(sp);
    for m in range(sp[0]):
      PafterOP[m]=np.matmul(OP,position[m,0:3].reshape((3,1))).reshape((1,3));
    print(OP)
    dp=analyzePH.diffposition(PafterOP,position,axis,sp[0]);
    if np.linalg.norm(dp)<1e-7:
      print("Check successful:----------------")
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
def findexactmatch(OP,axis,masslist,position):
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
    [PafterOPandT,Transform,PafterOP]=findexactmatch(OP,axis,masslist,position);
    transformlist.append(Transform);
  return transformlist;
if __name__=="__main__":
  pass;
else:
  pass;
