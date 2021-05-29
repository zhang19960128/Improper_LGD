import numpy as np
import symmop
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
    vOP=symmop.symoperate(symop,Tlist,Ev[gp[groupindex][0]],axis);
    re=findmodematchsubspace(vOP,Ev,gp[groupindex]);
    for j in range(20):
      vOPafter=symmop.symoperate(symop,Tlist,re[1],axis);
      re=findmodematchsubspace(vOPafter,Ev,gp[groupindex]);
    match=re[1].dot(vOPafter);
    if np.abs(match) > 0.80:
      matchcoeff.append(match);
      matchmode.append(re[1]);
    else:
      pass;
  return([matchcoeff,matchmode])
