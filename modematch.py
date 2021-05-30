import numpy as np
import symmop
def modematch(symop,Tlist,Ev,FREQ,axis,debug=0):
  matchcoeff=np.zeros(len(Ev));
  grouplist=groupmode(FREQ);
  for i in range(len(Ev)):
#    vOP=symmop.symoperate(symop,Tlist,Ev[i],axis);
    vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[i]);
    coeff=vOP.dot(Ev[i]);
    if np.abs(coeff) > 0.95:
      matchcoeff[i]=int(coeff/np.abs(coeff));
    else:
      groupid=findgroup(i,grouplist);
      if debug:
        print("group is: ",grouplist[groupid])
      Dimension=len(grouplist[groupid]);
      transmatrix=np.zeros((Dimension,Dimension));
      for j in range(Dimension):
        vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[grouplist[groupid][j]]);
        [coeff,matchpercentage]=findmodematchsubspace(vOP,Ev,grouplist[groupid],debug);
        if debug:
          print(symop)
        if np.abs(np.linalg.norm(coeff)-1.0) > 1e-1:
          pass
          if debug:
            print("It is not appropriately normed!!!")
        for k in range(Dimension):
          transmatrix[k][j]=coeff[k];
      matchcoeff[i]=0;
      if debug:
        print(transmatrix)
      for i in range(Dimension):
        matchcoeff[i]=matchcoeff[i]+transmatrix[i][i];
      if np.abs(matchpercentage) < 0.95:
        matchcoeff[i]=np.nan;
  return matchcoeff;
def findgroup(k,gplist):
  for i in range(len(gplist)):
    for j in range(len(gplist[i])):
      if np.abs(k - gplist[i][j]) < 1e-6:
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
def findmodematchsubspace(v1,v,subgroup,debug=0):
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
  if np.abs(subv.dot(v1)) < 0.95 and debug:
    print("CANNOT FIND SUBSPACE MATCH AFTER OPERATIONS!!!!!= matchcoeff:=",np.abs(subv.dot(v1)))
  return([coeff,subv.dot(v1)])
def modecheck(v):
  re=1;
  print("DM is:",len(v))
  for i in range(len(v)):
    for j in range(len(v)):
      if i==j:
        if np.abs(v[i].dot(v[j])-1.0) < 1e-6:
          re=re*1;
        else:
          re=re*0;
      else:
        if np.abs(v[i].dot(v[j])) < 1e-6:
          re=re*1;
        else:
          re=re*0;
  if np.abs(re-1) < 1e-6:
    print("Mode Orthonormal Successful");
  else:
    print("Mode Orthonormal Unsuccessful");
'''
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
def printmodetoscreen(Ev):
  for i in range(int(len(Ev)/3)):
    print("{0:10.7f} {1:10.7f} {2:10.7f}".format(Ev[3*i+0],Ev[3*i+1],Ev[3*i+2]));
'''
