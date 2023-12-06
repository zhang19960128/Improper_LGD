import numpy as np
import matplotlib.pyplot as plt
def getenergy(filename):
  fhandle=open(filename,'r');
  lines=fhandle.readlines();
  for i in range(len(lines)):
    if lines[i].find("ETOT")!=-1:
      ETOT=float(lines[i].split()[2]);
  return(ETOT)
xgrid=[];
ygrid=[];
length=20;
ygrid=np.zeros(length)
for i in range(20):
  xgrid.append(i);
  etot=getenergy("./T"+str(i)+"/PATH"+str(i)+".abo");
  ygrid[i]=etot;
  print('i=',i,ygrid[i])
plt.plot(xgrid,(ygrid-np.min(ygrid))*27.2,linewidth=4, markersize=12,label="6 modes");
plt.ylabel("Energy/(eV)")
plt.xlabel("PATH")
plt.xlim([0,20])
plt.show()
