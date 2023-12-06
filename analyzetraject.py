import FIXABI
import numpy as np
import matplotlib.pyplot as plt
def obtainH(filename):
  filehandle=open(filename,'r');
  lines=filehandle.readlines();
  length=len(lines);
  etot=0.0;
  for i in range(length):
    if lines[i].find("ETOT")!=-1:
      line=lines[i].split();
      etot=float(line[2]);
  return etot;
def obtainP():
  abi=FIXABI.abiout("./Pprofile/scf0","./Pprofile/scf0.abo","./Pprofile/dfpt0","./Pprofile/dfpt0.abo");
  scfinlist=[];
  scfoutlist=[];
  scfzeroinlist=[];
  scfzerooutlist=[];
  volume=np.linalg.det(abi.axis)*10**(-30);
  efieldabi=514220624373.482;
  echarge=1.602176634*10**(-19);
  HatoEv=27.2114;
  length=20;
  for i in range(length):
    scfinlist.append("./Pprofile/scf{0:d}".format(i));
    scfoutlist.append("./Pprofile/scf{0:d}.abo".format(i));
    scfzeroinlist.append("./Pprofile/itenoe{0:d}".format(i));
    scfzerooutlist.append("./Pprofile/itenoe{0:d}.abo".format(i));
  H0list=[];
  Hlist=[];
  polarlist=[];
  polarnoe=[];
  H0noe=[];
  efieldall=[];
  for i in range(length):
    polar=abi.obtainpolarization(scfoutlist[i]);
    if i==0:
      polarrefp=polar;
    polarref=abi.obtaindipolediffperiodtwo(scfinlist[0],scfoutlist[0],scfzerooutlist[0],scfinlist[i],scfoutlist[i],scfzerooutlist[i]);
    print("Polar reference is: ",polarref)
    polarlist.append((polarref+polarrefp)*57.21476);
    energy=obtainH(scfoutlist[i]);
    efield=abi.obtainefield(scfinlist[i]);
    efieldall.append(efield);
    H0=energy*HatoEv+np.dot(efield*efieldabi,(polar)*57.21476/echarge*volume);
    H0list.append(H0);
    Hlist.append(energy*HatoEv);
    polarref=abi.obtaindipolediffperiodtwo(scfinlist[0],scfoutlist[0],scfzerooutlist[0],scfzeroinlist[i],scfzerooutlist[i],scfzerooutlist[i]);
    polarnoe.append((polarref+polarrefp)*57.21476);
    energy=obtainH(scfzerooutlist[i]);
    H0noe.append(energy*HatoEv);
  pxlist=[];
  pylist=[];
  pzlist=[];
  for i in range(length):
    pxlist.append(polarlist[i][0]);
    pylist.append(polarlist[i][1]);
    pzlist.append(polarlist[i][2]);
    print(str(i)+" {0:12.5f} {1:12.5f} {2:12.5f} {3:12.5f} {4:12.5f} {5:12.5f} {6:12.5f} {7:12.5f}".format(polarlist[i][0],polarlist[i][1],polarlist[i][2],H0list[i],polarnoe[i][0],polarnoe[i][1],polarnoe[i][2],H0noe[i]))
  plt.rc('font', size=16)
  plt.plot(pylist,(H0list-min(H0list))*1000,'ro-',linewidth=4, markersize=12,label= "m2")
  plt.xlabel("Polarization(C/m^2)")
  plt.ylabel("Energy(mev)")
  plt.show()
  return([pxlist,pylist,pzlist])
