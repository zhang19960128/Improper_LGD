import numpy as np
import symmop
def modematch(symop,Tlist,Ev,FREQ,axis):
  matchcoeff=np.zeros(len(Ev),dtype=int);
  grouplist=groupmode(FREQ);
  for i in range(len(Ev)):
#    vOP=symmop.symoperate(symop,Tlist,Ev[i],axis);
    vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[i]);
    coeff=vOP.dot(Ev[i]);
    if np.abs(coeff) > 0.95:
      matchcoeff[i]=int(coeff/np.abs(coeff));
    else:
      groupid=findgroup(i,grouplist);
      Dimension=len(grouplist[groupid]);
      transmatrix=np.zeros((Dimension,Dimension));
      for j in range(Dimension):
        vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[grouplist[groupid][j]]);
      matchcoeff[i]=0;
  return matchcoeff;
def findgroup(k,gplist):
  for i in range(len(gplist)):
    for j in range(len(gplist[i])):
      if np.abs(k - gplist[i][j]) < 1e-6;
        return i;
  return -1;
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
  if np.abs(subv.dot(v1)) < 0.90:
    print("CANNOT FIND SUBSPACE MATCH AFTER OPERATIONS!!!!!")
  return(coeff)
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
    if np.abs(match) > 0.90:
      matchcoeff.append(match);
      matchmode.append(re[1]);
    else:
      pass;
  return([matchcoeff,matchmode])
