import numpy as np
import inputfiles
import sequence
def readxsf(filename,natom):
  f=open(filename,'r');
  lines=f.readlines();
  position=np.zeros((natom,3));
  for i in range(8,len(lines)):
    line=lines[i].split();
    for j in range(3):
      position[i-8][j]=float(line[j]);
  return position
#Rposition-BasePosition#
def diffP(Rposition,BasePosition):
  diffp=np.zeros(np.shape(BasePosition))
  for i in range(len(BasePosition)):
    minimum=1e8;
    for j in range(len(Rposition)):
      dp=Rposition[j]-BasePosition[i]-np.round(Rposition[j]-BasePosition[i]);
      if np.linalg.norm(dp) < minimum:
        diffp[i]=np.copy(dp);
        minimum=np.linalg.norm(dp);
  return diffp
def projection(eigenvector,dp,masslist):
  dptoeig=np.copy(dp);
  for i in range(len(dptoeig)):
    dptoeig[i]=np.sqrt(masslist[i])*dp[i];
  sp=np.shape(dptoeig);
  dp1D=dp.reshape(sp[0]*sp[1]);
  Dimension=len(eigenvector);
  coeff=np.zeros(Dimension);
  for i in range(Dimension):
    coeff[i]=eigenvector[i].dot(dp1D);
  return coeff;
def analyzetransition(vGa,vM):
  Ophase=readxsf(inputfiles.EXPORTHOPOSCAR,inputfiles.natomExp);
  Tphase=readxsf(inputfiles.EXPTETRAPOSCAR,inputfiles.natomExp);
  Cphase=readxsf(inputfiles.EXPCUBICPOSCAR,inputfiles.natomExp);
  [masslist,namelist]=sequence.sequence(inputfiles.natomExp)
  dpTC=diffP(Tphase,Cphase);
  dpOC=diffP(Ophase,Cphase);
  pTCtoGa=projection(vGa,dpTC,masslist);
  pTCtoM=projection(vM,dpTC,masslist);
  pOCtoGa=projection(vGa,dpOC,masslist);
  pOCtoM=projection(vM,dpOC,masslist);
  return([pTCtoGa,pTCtoM,pOCtoGa,pOCtoM])
