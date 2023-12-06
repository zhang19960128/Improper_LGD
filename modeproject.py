import numpy as np
import positdiff
def modeproject(file1,file2,natom,axis,masslist,vmerge):
  dp=positdiff.dp(natom,axis,file1,file2);
  dp1D=dp.reshape(3*natom);
  modepro=np.zeros(3*natom);
  mass3D=[];
  eigendisplist=[];
  for i in range(natom):
    for j in range(3):
      mass3D.append(masslist[i]);
  for i in range(3*natom):
    eigendisp=np.copy(vmerge[i]);
#    for j in range(3*natom):
#      eigendisp[j]=eigendisp[j]/np.sqrt(mass3D[j]);
#    eigendisp=eigendisp/np.linalg.norm(eigendisp);
    modepro[i]=(dp1D*(np.sqrt(np.array(mass3D)))).dot(vmerge[i]);
#    modepro[i]=dp1D.dot(eigendisp);
    eigendisplist.append(eigendisp);
#    nm=np.linalg.norm(vmerge[i][12:15]);
#    if nm > 1e-5:
#      modepro[i]=modepro[i]/nm;
#    else:
#      modepro[i]=modepro[i]/1.0;
  return [modepro,eigendisplist];
'''
  for i in range(3*natom):
    print("Mode ID is:=",i,"{0:10.7f}".format((dp1D*(np.sqrt(np.array(mass3D)))).dot(vmerge[i])))
'''

