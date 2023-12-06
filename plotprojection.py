import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import modeproject
import positdiff
import analyzetraject
def projectpath(natom,axis,masslist,vmerge):
  m0=[];
  m1=[];
  m2=[];
  m3=[];
  m4=[];
  m5=[];
  m6=[];
  m7=[];
  m8=[];
  m9=[];
  m10=[];
  m11=[];
  x=[];
  y=[];
  z=[];
  m15=[];
  m16=[];
  m17=[];
  m18=[];
  m19=[];
  m20=[];
  m21=[];
  m22=[];
  m23=[];
  m24=[];  
  m25=[];
  m26=[];
  m27=[];
  m28=[];
  m29=[];
  length=20;
  for i in range(length):
    [modepro,eigendisplist]=modeproject.modeproject("./ABIOUTPUT/dfpt0666","./Pprofile/scf"+str(i),natom,axis,masslist,vmerge);
    m0.append(modepro[0]);
    m1.append(modepro[1]);
    m2.append(modepro[2]);
    m3.append(modepro[3]);
    m4.append(modepro[4]);
    m5.append(modepro[5]);
    m6.append(modepro[6]);
    m7.append(modepro[7]);
    m8.append(modepro[8]);
    m9.append(modepro[9]);
    m10.append(modepro[10]);
    m11.append(modepro[11]);
    x.append(modepro[12]);
    y.append(modepro[13]);
    z.append(modepro[14]);
    m15.append(modepro[15]);
    m16.append(modepro[16]);
    m17.append(modepro[17]);
    m18.append(modepro[18]);
    m19.append(modepro[19]);
    m20.append(modepro[20]);
    m21.append(modepro[21]);
    m22.append(modepro[22]);
    m23.append(modepro[23]);
    m24.append(modepro[24]);
    m25.append(modepro[25]);
    m26.append(modepro[26]);
    m27.append(modepro[27]);
    m28.append(modepro[28]);
    m29.append(modepro[29]);
  [px,py,pz]=analyzetraject.obtainP();
  basep=positdiff.readposition("./Pprofile/dfpt0666",axis,natom);
  for i in range(length):
    pchange=m0[i]*eigendisplist[0];
    pchange=m1[i]*eigendisplist[1]+pchange;
    pchange=m2[i]*eigendisplist[2]+pchange;
    pchange=m3[i]*eigendisplist[3]+pchange;
    pchange=m4[i]*eigendisplist[4]+pchange;
    pchange=m5[i]*eigendisplist[5]+pchange;
    pchange=m6[i]*eigendisplist[6]+pchange;
    pchange=m7[i]*eigendisplist[7]+pchange;
    pchange=m8[i]*eigendisplist[8]+pchange;
    pchange=m9[i]*eigendisplist[9]+pchange;
    pchange=m10[i]*eigendisplist[10]+pchange;
    pchange=m11[i]*eigendisplist[11]+pchange;
    pchange=x[i]*eigendisplist[12]+pchange;
    pchange=y[i]*eigendisplist[13]+pchange;
    pchange=z[i]*eigendisplist[14]+pchange;
    pchange=m15[i]*eigendisplist[15]+pchange;
    pchange=m16[i]*eigendisplist[16]+pchange;
    pchange=m17[i]*eigendisplist[17]+pchange;
    pchange=m18[i]*eigendisplist[18]+pchange;
    pchange=m19[i]*eigendisplist[19]+pchange;
    pchange=m20[i]*eigendisplist[20]+pchange;
    pchange=m21[i]*eigendisplist[21]+pchange;
    pchange=m22[i]*eigendisplist[22]+pchange;
    pchange=m23[i]*eigendisplist[23]+pchange;
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
  plt.rc('font', size=16)   
  plt.plot(py,m2,"-bo",linewidth=4, markersize=12, label=r'm$_{02}$')  
  plt.plot(py,m11,"-yo",linewidth=4, markersize=12,label=r'm$_{25}$')
  plt.plot(py,y,"-ro",linewidth=4, markersize=12,  label=r'$y$')
  plt.plot(py,m15,"-ko",linewidth=4, markersize=12,label=r'm$_{40}$')
  plt.plot(py,m19,"-bD",linewidth=4, markersize=12,label=r'm$_{51}$')
  plt.plot(py,m21,"-gD",linewidth=4, markersize=12,label=r'm$_{53}$')
  plt.plot(py,m22,"-rD",linewidth=4, markersize=12,label=r'm$_{54}$')
  plt.plot(py,m29,"-yD",linewidth=4, markersize=12,label=r'm$_{65}$')
  coefffile=open("COEFF.dat",'w');
  for i in range(len(m2)):
    coefffile.write("{0:10.7f} {1:10.7f} {2:10.7f} {3:10.7f} {4:10.7f} {5:10.7f} {6:10.7f} {7:10.7f}\n".format(m2[i],m11[i],y[i],m15[i],m19[i],m21[i],m22[i],m29[i]))
  plt.legend(loc="lower left");
  plt.xlabel("Polarization($C/m^2$)")
  plt.ylabel("Projection on Mode")
  plt.rc('axes',labelsize=100)   
  plt.rc('axes', titlesize=100) 
  plt.show();
