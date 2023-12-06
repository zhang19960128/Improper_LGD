import numpy as np
import matplotlib.pyplot as plt
import modeproject
import positdiff
import analyzetraject
def projectpath(natom,axis,masslist,vmerge):
  a=[];
  b=[];
  c=[];
  x=[];
  y=[];
  z=[];
  e=[];
  f=[];
  g=[];
  k=[];
  m=[];
  n=[];
  m3=[];
  m4=[];
  m5=[];
  m6=[];
  m7=[];
  m8=[];\
  m9=[];
  m10=[];
  m11=[];
  m15=[];
  m16=[];
  m17=[];
  m24=[];  
  m25=[];
  m26=[];
  m27=[];
  m28=[];
  m29=[];
  length=20;
  for i in range(length):
    [modepro,eigendisplist]=modeproject.modeproject("./ABIOUTPUT/dfpt0666","./Pprofile/scf"+str(i),natom,axis,masslist,vmerge);
    a.append(modepro[0]);
    b.append(modepro[1]);
    c.append(modepro[2]);
    x.append(modepro[12]);
    y.append(modepro[13]);
    z.append(modepro[14]);
    e.append(modepro[18]);
    f.append(modepro[19]);
    g.append(modepro[20]);
    k.append(modepro[21]);
    m.append(modepro[22]);
    n.append(modepro[23]);
    m3.append(modepro[3]);
    m4.append(modepro[4]);
    m5.append(modepro[5]);
    m6.append(modepro[6]);
    m7.append(modepro[7]);
    m8.append(modepro[8]);
    m9.append(modepro[9]);
    m10.append(modepro[10]);
    m11.append(modepro[11]);
    m15.append(modepro[15]);
    m16.append(modepro[16]);
    m17.append(modepro[17]);
    m24.append(modepro[24]);
    m25.append(modepro[25]);
    m26.append(modepro[26]);
    m27.append(modepro[27]);
    m28.append(modepro[28]);
    m29.append(modepro[29]);
  [px,py,pz]=analyzetraject.obtainP();
  basep=positdiff.readposition("./Pprofile/dfpt0666",axis,natom);
  for i in range(length):
    pchange=a[i]*eigendisplist[0];
    pchange=b[i]*eigendisplist[1]+pchange;
    pchange=c[i]*eigendisplist[2]+pchange;
#    pchange=x[i]*eigendisplist[12]+pchange;
#    pchange=y[i]*eigendisplist[13]+pchange;
#    pchange=z[i]*eigendisplist[14]+pchange;
    pchange=e[i]*eigendisplist[18]+pchange;
    pchange=f[i]*eigendisplist[19]+pchange;
    pchange=g[i]*eigendisplist[20]+pchange;
    pchange=k[i]*eigendisplist[21]+pchange;
    pchange=m[i]*eigendisplist[22]+pchange;
    pchange=n[i]*eigendisplist[23]+pchange;
    pchange=m3[i]*eigendisplist[3]+pchange;
    pchange=m4[i]*eigendisplist[4]+pchange;
    pchange=m5[i]*eigendisplist[5]+pchange;
    pchange=m6[i]*eigendisplist[6]+pchange;
    pchange=m7[i]*eigendisplist[7]+pchange;
    pchange=m8[i]*eigendisplist[8]+pchange;
    pchange=m9[i]*eigendisplist[9]+pchange;
    pchange=m10[i]*eigendisplist[10]+pchange;
    pchange=m11[i]*eigendisplist[11]+pchange;
    pchange=m15[i]*eigendisplist[15]+pchange;
    pchange=m16[i]*eigendisplist[16]+pchange;
    pchange=m17[i]*eigendisplist[17]+pchange;
    pchange=m24[i]*eigendisplist[24]+pchange;
    pchange=m25[i]*eigendisplist[25]+pchange;
    pchange=m26[i]*eigendisplist[26]+pchange;
    pchange=m27[i]*eigendisplist[27]+pchange;
    pchange=m28[i]*eigendisplist[28]+pchange;
    pchange=m29[i]*eigendisplist[29]+pchange;
    preshape=pchange.reshape((natom,3))
    for j in range(natom):
      preshape[j]=preshape[j]/np.sqrt(masslist[j]);
    newp=basep+preshape;
    positdiff.writenewscfnoe("./Pprofile/scf0template",natom,axis,newp,"/workspace/jiahaoz/HZO/LG_Model/GroundP/HfO4/WaveVector/process/Improper_LGD/PATH/PATH"+str(i))
'''
  plt.rc('font', size=16)   
  plt.plot(py,a,"-bo",linewidth=4, markersize=12,label="AFE(C) kx Order Para=a")
  plt.plot(py,b,"-yo",linewidth=4, markersize=12,label="AFE(C) ky Order Para=b")
  plt.plot(py,c,"-ro",linewidth=4, markersize=12,label="AFE(C) kz Order Para=c")
  plt.plot(py,y,"-ko",linewidth=4, markersize=12,label="Fey   Order Para=y")
  plt.plot(py,e,"-bD",linewidth=4, markersize=12,label="Kx(z) Order Para=e")
  plt.plot(py,f,"-gD",linewidth=4, markersize=12,label="Kx(y) Order Para=f")
  plt.plot(py,g,"-rD",linewidth=4, markersize=12,label="Ky(x) Order Para=g")
  plt.plot(py,k,"-yD",linewidth=4, markersize=12,label="Ky(z) Order Para=k")
  plt.plot(py,m,"-kD",linewidth=4, markersize=12,label="Kz(x) Order Para=m")
  plt.plot(py,n,"-mD",linewidth=4, markersize=12,label="Kz(y)  Order Para=n")
  print(f[-1])
  plt.legend(loc="lower left");
  plt.xlabel("Polarization(C/m^2)")
  plt.ylabel("Projection on Mode")
  plt.rc('axes',labelsize=100)   
  plt.rc('axes', titlesize=100) 
  plt.show();
'''
