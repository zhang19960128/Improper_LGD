import numpy as np
import math
import matplotlib.pyplot as plt
def writenewscfnoe(scfin,natoms,axis,atomposition,filename):
  scffiles=open(scfin,'r');
  newfilename=open(filename,'w');
  lines=scffiles.readlines();
  skiplines=0;
  for i in range(len(lines)):
    if lines[i].find('efield')!=-1:
      skiplines=skiplines+1;
    elif lines[i].find('berryopt')!=-1:
      skiplines=skiplines+1;
    elif lines[i].find("xred")!=-1:
      newfilename.write("xred ");
      for j in range(natoms):
        temp="";
        atomp=np.matmul(atomposition[j][0:3],np.linalg.inv(axis.transpose()));
        for k in range(3):
          temp=temp+' {:12.8f}'.format(atomp[k]);
        if j==0:
          newfilename.write(temp+'\n');
        else:
          newfilename.write("    "+temp+'\n');
        skiplines=skiplines+1;
    elif lines[i].find("xcart")!=-1:
      if lines[i+natoms-1].lower().find('bohr')!=-1:
        autoA=0.529;
        tag='bohr';
      elif lines[i+natoms-1].lower().find('angstrom')!=-1:
        autoA=1.0;
        tag='angstrom';
      newfilename.write("xcart ");
      for j in range(natoms):
        temp="";
        atomp=atomposition[j]/autoA;
        for k in range(3):
          temp=temp+' {:12.8f}'.format(atomp[k]);
        if j==0:
          newfilename.write(temp+'\n');
        elif j==natoms-1:
          newfilename.write(temp+' '+tag+'\n');
        else:
          newfilename.write("    "+temp+'\n');
        skiplines=skiplines+1;
    if skiplines >0:
      skiplines=skiplines-1;
      continue;
    else:
      newfilename.write(lines[i]);
  newfilename.close();
  scffiles.close();
def readposition(scfin,axis,natoms):
    scfinput=open(scfin,'r');
    lines=scfinput.readlines();
    length=len(lines);
    atomp=np.zeros([natoms,3]);
    atomreturn=np.zeros([natoms,3]);
    for t in range(length):
      if lines[t].find("xred")!=-1:
        for i in range(natoms):
          line=lines[t+i].split();
          for j in range(-3,0,1):
            atomp[i][j]=float(line[j]);
        for i in range(natoms):
          for j in range(3):
            atomreturn[i]=atomreturn[i]+axis[j]*atomp[i][j];
      if lines[t].find("xcart")!=-1:
        for i in range(natoms):
          line=lines[t+i].split();
          if len(lines[t+natoms-1].split())==3:
            autoA=0.52917721067121;
            for i in range(natoms):
              line=lines[t+i].split();
              for j in range(-3,0,1):
                atomreturn[i][j]=float(line[j])*autoA;
          elif len(lines[t+natoms-1].split())==4:
            if lines[t+natoms-1].lower().find("bohr"):
              autoA=0.52917721067121;
            elif lines[t+natoms-1].lower().find("angstrom"):
              autoA=1.0;
            line=lines[t+0].split();
            for j in range(-3,0,1):
              atomreturn[0][j]=float(line[j])*autoA;
            for i in range(1,natoms):
              line=lines[t+i].split();
              for j in range(3):
                atomreturn[i][j]=float(line[j])*autoA;
    return atomreturn;
def writeinterpolate(natoms,axis,atomposition,filein,filename):
  files=open(filein,'r');
  lines=files.readlines();
  files.close();
  newfilename=open(filename,'w');
  skiplines=0;
  for i in range(len(lines)):
    if lines[i].find("xred")!=-1:
      newfilename.write("xred ");
      for j in range(natoms):
        temp="";
        atomp=np.matmul(atomposition[j][0:3],np.linalg.inv(axis.transpose()));
        for k in range(3):
          temp=temp+' {:12.8f}'.format(atomp[k]);
        if j==0:
          newfilename.write(temp+'\n');
        else:
          newfilename.write("    "+temp+'\n');
        skiplines=skiplines+1;
    elif lines[i].find("xcart")!=-1:
      if lines[i+natoms-1].lower().find('bohr')!=-1:
        autoA=0.529;
      elif lines[i+natoms-1].lower().find('angstrom')!=-1:
        autoA=1.0;
      newfilename.write("xcart ");
      for j in range(natoms):
        temp="";
        atomp=atomposition[j]/autoA;
        for k in range(3):
          temp=temp+' {:12.8f}'.format(atomp[k]);
        if j==0:
          newfilename.write(temp+'\n');
        else:
          newfilename.write("    "+temp+'\n');
        skiplines=skiplines+1;
    if skiplines >0:
      skiplines=skiplines-1;
      continue;
    else:
      newfilename.write(lines[i]);
  newfilename.close();
def distance(p1,p2,axis):
  dis=0.0;
  for i in range(3):
    temp=p1[i]-p2[i]-round((p1[i]-p2[i])/axis[i][i])*axis[i][i];
    dis=dis+temp*temp;
  return math.sqrt(dis)
def obtainpolarization(scfoutfile):
  scfout=open(scfoutfile,'r');
  lines=scfout.readlines();
  polar=np.zeros(3);
  for t in range(len(lines)):
    if lines[t].find("Polarization in cartesian coordinates (a.u.)")!=-1:
      line=lines[t+5].split();
      for j in range(-3,0,1):
        polar[j]=float(line[j]);
  scfout.close();
  return polar; # the units of Polarization is e/bohr^2;
def dp(natom,axis,file1,file2):
  psymm=readposition(file1,axis,natom);
  pground=readposition(file2,axis,natom);
  for i in range(1,natom):
    for j in range(3):
      psymm[i][j]=psymm[i][j]-psymm[0][j];
      pground[i][j]=pground[i][j]-pground[0][j];
  for j in range(3):
    psymm[0][j]=0.0;
    pground[0][j]=0.0;
  pgroundaftersort=np.zeros((natom,3));
  dist=np.zeros(natom);
  for i in range(0,natom):
    for j in range(0,natom):
      dist[j]=distance(psymm[i],pground[j],axis);
    sortlist=np.argsort(dist);
    for k in range(3):
      pgroundaftersort[i][k]=pground[sortlist[0]][k];
      pgroundaftersort[i][k]=pgroundaftersort[i][k]-round((pgroundaftersort[i][k]-psymm[i][k])/axis[k][k])*axis[k][k];
  return (pgroundaftersort-psymm);
