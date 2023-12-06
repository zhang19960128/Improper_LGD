import numpy as np
import IOmode
import symmop
import check
def modematch(symop,Tlist,Ev,FREQ,axis):
  matchcoeff=np.zeros(len(Ev));
  grouplist=groupmode(FREQ);
  for i in range(len(Ev)):
    vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[i]);#vOP=symmop.symoperate(symop,Tlist,Ev[i],axis);
    coeff=vOP.dot(Ev[i]);
    if np.abs(coeff) > 0.95:
      matchcoeff[i]=coeff/np.abs(coeff);
    else:
      groupid=findgroup(i,grouplist);
      Dimension=len(grouplist[groupid]);
      transmatrix=np.zeros((Dimension,Dimension));
      for j in range(Dimension):
        vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[grouplist[groupid][j]]);
        [coeff,matchpercentage]=findmodematchsubspace(vOP,Ev,grouplist[groupid]);
        if np.abs(np.linalg.norm(coeff)-1.0) > 1e-1:
          print("NOT Normalized!",np.linalg.norm(coeff),'coeff:=',coeff,'matchpercentage is:',matchpercentage,'groupid is:=',groupid);
        else:
          pass
        for k in range(Dimension):
          transmatrix[k][j]=coeff[k];
        if np.abs(matchpercentage) < 0.95:
          transmatrix[j][j]=np.nan;
      matchcoeff[i]=np.trace(transmatrix);
      if np.abs(matchpercentage) < 0.95:
        matchcoeff[i]=np.nan;
  return matchcoeff;
def modematchsecond(symop,Tlist,Ev,FREQ,axis):
  matchcoeff=np.zeros(len(Ev));
  grouplist=groupmode(FREQ);
  for i in range(len(Ev)):
    groupid=findgroup(i,grouplist);
    if len(grouplist[groupid])==1:
      vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[i]);#vOP=symmop.symoperate(symop,Tlist,Ev[i],axis);
      coeff=vOP.dot(Ev[i]);
      if np.abs( np.abs(coeff) - 1 ) < 1e-3:
        matchcoeff[i]=coeff;
      else:
        print("FATAL ERROR, vector is not mapped in subspace");
    else:
      groupid=findgroup(i,grouplist);
      Dimension=len(grouplist[groupid]);
      transmatrix=np.zeros((Dimension,Dimension));
      for m in range(Dimension):
        for n in range(Dimension):
          vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[grouplist[groupid][n]]);
          transmatrix[m][n]=vOP.dot(Ev[grouplist[groupid][m]]);
      if( np.linalg.norm(transmatrix-np.identity(Dimension)) < 1e-6 ):
        matchcoeff[i]=1.0;
        matchcoeff[i]=np.trace(transmatrix);
      elif( np.linalg.norm(transmatrix+np.identity(Dimension)) < 1e-6 ):
        matchcoeff[i]=-1.0;
        matchcoeff[i]=np.trace(transmatrix);
      else:
        matchcoeff[i]=np.trace(transmatrix);
      for m in range(Dimension):
        coeff=transmatrix[:,m].reshape(Dimension);
        add=np.zeros(len(Ev[grouplist[groupid][m]]));
        for n in range(Dimension):
          add=add+coeff[n]*Ev[grouplist[groupid][n]];
        vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[grouplist[groupid][m]]);
        difference=np.linalg.norm(vOP-add);
        if difference > 1e-5:
          print("FATAL ERROR, not in subspace");
        else:
          pass;
  return matchcoeff;
def findgroup(k,gplist):
  for i in range(len(gplist)):
    for j in range(len(gplist[i])):
      if np.abs(k - gplist[i][j]) < 1e-6:
        return i;
  print("!!!!!!!!!!!!!!!!!! FATAL ERROR CANNOT FIND GROUP !!!!!!!!!!!!!!!")
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
      if np.abs(w[groupset[i][0]]-w[modeset[0]])/np.abs(w[groupset[i][0]]) < 5*1e-3:
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
  if np.abs(subv.dot(v1)) < 0.95:
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
          print("Self Not orthogonal")
      else:
        if np.abs(v[i].dot(v[j])) < 1e-6:
          re=re*1;
        else:
          re=re*0;
          print("Undiag is not orthogonal: ",v[i].dot(v[j]),i,j)
  if np.abs(re-1) < 1e-4:
    print("Mode Orthonormal Successful");
  else:
    print("Mode Orthonormal Unsuccessful");
    for i in range(len(v)):
      for j in range(len(v)):
        if i!=j:
          print("UnDiag=({0:2d},{1:2d} {2:10.7f})".format(i,j,v[i].dot(v[j])));
        else:
          print("Diag=({0:2d},{1:2d} {2:10.7f})".format(i,j,v[i].dot(v[j])));
def matchire(masslist,axis,position,modeid,mode,FREQ,opmatrix,opmatrixinrep):
  grouplist=groupmode(FREQ);
  groupindex=findgroup(modeid,grouplist);
  Tlist=check.check(opmatrix,axis,masslist,position);
  Opsubspace=[];
  for i in range(len(opmatrix)):
    vOP=symmop.symoperate_matrix_op(opmatrix[i],Tlista[i],mode[modeid]);
    if np.abs(np.abs(vOP.dot(mode[modeid]))-1.0)<1e-6:
      Opsubspace.append(vOP.dot(mode[modeid]));
    else:
      groupid=findgroup(i,grouplist);
      Dimension=len(grouplist[groupid]);
      transmatrix=np.zeros((Dimension,Dimension));
      for m in range(Dimension):
        for n in range(Dimension):
          vOP=symmop.symoperate_matrix_op(symop,Tlist,Ev[grouplist[groupid][n]]);
          transmatrix[m][n]=vOP.dot(Ev[grouplist[groupid][m]]);
