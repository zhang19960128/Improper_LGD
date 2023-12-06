import numpy as np
import cmath
Ha=2*13.59*1.602176634*10**(-19);
bohr=0.529177*10**-10;
aumass=1.660539066*10**-27;
def printmode(Hposition,namelist,masslist,axis,v,modename):
  natom=len(namelist);
  for i in range(len(v)):
    f=open(modename+str(i)+".xsf",'w');
    f2=open(modename+str(i),'w');
    for j in range(3*natom):
      f2.write("{0:10.7f} ".format(v[i][j]));
    f2.close();
    f.write("PRIMVEC\n");
    for j in range(3):
      for k in range(3):
        f.write(str(axis[j][k])+" ");
      f.write("\n");
    f.write("PRIMCOORD\n");
    f.write(str(natom)+" 1\n");
    for j in range(natom):
      f.write(namelist[j]+" {0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f} {4:10.7f} {5:10.7f}\n".format(Hposition[j][0],Hposition[j][1],Hposition[j][2],v[i][3*j+0]/5.0/np.sqrt(masslist[j]),v[i][3*j+1]/5.0/np.sqrt(masslist[j]),v[i][3*j+2]/5.0/np.sqrt(masslist[j])))
    f.close();
def printmodefreq(w,natom):
  for i in range(3*natom):
    print(cmath.sqrt(w[i]*Ha/bohr/bohr/aumass)*(10**-12)*33.35641/2/3.141592653)
def printvector(vec):
  temp="";
  length=len(vec);
  for i in range(length):
    temp=temp+" "+"{0:7.3f}".format(vec[i]);
  print(temp)
def printMatrix(s):
  print('-----------------------------------------')
  for i in range(len(s)):
    for j in range(len(s[0])):
      print("{0:4.3f} ".format(s[i][j]),end="")
    print('\n')
  print('-----------------------------------------')
