import numpy as np
import symmop
import modematch
def irrep(symoplist,Ev,FREQ,axis):
  grouplist=modematch.groupmode(FREQ);
  for groupid in range(len(grouplist)):
    groupirrepfile=open("IREP"+str(groupid)+'.dat','w');
    for j in range(len(symoplist)):
      if len(grouplist[groupid])==1:
        symop=symoplist[j];
        vOP=symmop.OPdotV(symop,Ev[grouplist[groupid][0]]);#vOP=symmop.symoperate(symop,Tlist,Ev[i],axis);
        coeff=vOP.dot(Ev[0]);
        if np.abs( np.abs(coeff) - 1 ) < 1e-3:
          pass;
        else:
          print("FATAL ERROR, vector is not mapped in subspace");
      else:
        symop=symoplist[j];
        Dimension=len(grouplist[groupid]);
        transmatrix=np.zeros((Dimension,Dimension));
        for m in range(Dimension):
          for n in range(Dimension):
            vOP=symmop.OPdotV(symop,Ev[grouplist[groupid][n]]);
            transmatrix[m][n]=vOP.dot(Ev[grouplist[groupid][m]]);
        for m in range(Dimension):
          coeff=transmatrix[:,m].reshape(Dimension);
          add=np.zeros(len(Ev[grouplist[groupid][m]]));
          for n in range(Dimension):
            add=add+coeff[n]*Ev[grouplist[groupid][n]];
          vOP=symmop.OPdotV(symop,Ev[grouplist[groupid][m]]);
          difference=np.linalg.norm(vOP-add);
          if difference > 1e-5:
            print("FATAL ERROR, not in subspace","OP id=",j,"Group ID:=",groupid);
          else:
            pass;
      for m in range(np.shape(transmatrix)[0]):
        for n in range(np.shape(transmatrix)[1]): 
          groupirrepfile.write("{0:10.7f} ".format(transmatrix[m][n]));
        groupirrepfile.write("\n");
def mergeirrep(symopnum,FREQ,mergelist):
  grouplist=modematch.groupmode(FREQ);
  TotalDm=0;
  Dmlist=[];
  for i in range(len(mergelist)):
    TotalDm=TotalDm+len(grouplist[mergelist[i]]);
    Dmlist.append(len(grouplist[mergelist[i]]));
  Ireplist=[];
  for i in range(len(mergelist)):
    irepdat=np.loadtxt("IREP"+str(mergelist[i])+'.dat')
    Ireplist.append(irepdat);
  totalfilestr="";
  for i in range(len(mergelist)):
    totalfilestr=totalfilestr+"+"+str(mergelist[i]);
  mergeirepfhandle=open("IREP"+totalfilestr+'.dat','w');
  totalmatrix=[];
  for i in range(symopnum):
    matrixforprint=np.zeros((TotalDm,TotalDm));
    for m in range(len(mergelist)):
      sumdim=0;
      for n in range(m):
        sumdim=Dmlist[n]+sumdim;
      matrixforprint[sumdim:(sumdim+Dmlist[m]),sumdim:(sumdim+Dmlist[m])]=np.copy(Ireplist[m][i*Dmlist[m]:(i+1)*Dmlist[m],:]);
    totalmatrix.append(matrixforprint);
#  total=[];
#  totalset=[];
#  for i in range(len(totalmatrix)):
#    for j in range(len(totalmatrix)):
#      total.append(np.matmul(totalmatrix[i],totalmatrix[j]));
#  for i in range(len(total)):
#    totalset.append(tuple(total[i].reshape(TotalDm*TotalDm)));
#  totalaftermerge=set(totalset);
#  totalaftermerge=list(totalaftermerge);
#  groupmatrix=[];
#  for i in range(len(totalaftermerge)):
#    element=np.array(totalaftermerge[i]);
#    matrix=element.reshape((TotalDm,TotalDm));
  for i in range(len(totalmatrix)):
    for m in range(TotalDm):
      for n in range(TotalDm):
        mergeirepfhandle.write("{0:10.4f} ".format(totalmatrix[i][m][n]));
      mergeirepfhandle.write("\n");
  mergeirepfhandle.close();
