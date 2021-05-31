import numpy as np
import math
Ha=2*13.59*1.602176634*10**(-19);
bohr=0.529177*10**-10;
aumass=1.660539066*10**-27;
emass=9.1093837015*10**(-31);
def obtaindyn(natoms,dynfile):
   dynmatrix=np.zeros((natoms,natoms,3,3));
   dyn=open(dynfile,'r');
   lines=dyn.readlines(); # note that the dynamical matrix obtain from abinit calculation is acutally force constant matrix Units Ha/bohr/bohr
   for t in range(len(lines)):
     if lines[t].find("Dynamical matrix, in cartesian coordinates")!=-1:
       tick=0;
       for i in range(natoms):
         for m in range(3):
           for j in range(natoms):
             for n in range(3):
               line=lines[t+tick+5].split();
               for s in range(4):
                 line[s]=int(line[s]);
               line[4]=float(line[4]);
               dynmatrix[line[1]-1][line[3]-1][line[0]-1][line[2]-1]=line[4];
               tick=tick+1; 
           tick=tick+1; # remove one blank
   dyn.close();
   return dynmatrix;
def readaxis(filename):
  files=open(filename,'r');
  lines=files.readlines();
  axis=np.zeros((3,3));
  scale=np.zeros(3);
  for i in range(len(lines)):
    if lines[i].find("acell")!=-1:
      if lines[i].lower().find("bohr")!=-1:
        autoA=0.529;
      elif lines[i].lower().find("angstrom")!=-1:
        autoA=1.0;
      else:
        autoA=0.529;
      line=lines[i].split();
      for j in range(3):
        scale[j]=float(line[j+1])*autoA;
    if lines[i].find("rprim")!=-1:
      for j in range(3):
        line=lines[i+j].split();
        for k in range(-3,0,1):
          axis[j][k]=float(line[k]);
  reaxis=np.zeros((3,3));
  for i in range(3):
    reaxis[i]=scale[i]*axis[i];
  return reaxis;
def readposition(axis,scfin,natoms):
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
def washback(position,natom):
  for i in range(natom-1,-1,-1):
    for j in range(3):
      position[i,j]=position[i,j]-position[0,j];
def distance(p1,p2,axis):
  dis=0.0;
  for i in range(3):
    temp=p1[i]-p2[i]-round((p1[i]-p2[i])/axis[i][i])*axis[i][i];
    dis=dis+temp*temp;
  return math.sqrt(dis)
def diffposition(PA,PB,axis,natom):
  dpaftersort=np.zeros((natom,3));
  PAaftersort=np.zeros((natom,3));
  dist=np.zeros(natom);
  for i in range(0,natom):
    for j in range(0,natom):
      dist[j]=distance(PB[i],PA[j],axis);
    sortlist=np.argsort(dist);
    for k in range(3):
      PAaftersort[i][k]=PA[sortlist[0]][k];
      temp=(PAaftersort[i][k]-PB[i][k])/axis[k][k];
      PAaftersort[i][k]=PB[i][k]+(temp-round(temp))*axis[k][k];
  return PAaftersort-PB;
def readsymmetry(filename):
  f=open(filename,'r');
  lines=f.readlines();
  length=len(lines);
  symmoperations=[];
  firstfind=1;
  for i in range(length):
    if lines[i].find("symrel")!=-1 and firstfind==1:
      k=0;
      symatrix=np.zeros((2,3,3));
      while len(lines[i+k].split())==18 or k==0 or len(lines[i+k].split())==9:
        listline=lines[i+k].split();
        if len(listline)==19 and k==0:
          for t in range(2):
            for m in range(3):
              for n in range(3):
                symatrix[t][m][n]=float(listline[t*3*3+3*m+n+1]);
            symmoperations.append(np.copy(symatrix[t,:,:]));
        if len(listline)==18:
          for t in range(2):
            for m in range(3):
              for n in range(3):
                symatrix[t][m][n]=float(listline[t*3*3+3*m+n]);
            symmoperations.append(np.copy(symatrix[t,:,:]));
        if len(listline)==10:
          for t in range(1):
            for m in range(3):
              for n in range(3):
                symatrix[t][m][n]=float(listline[t*3*3+3*m+n+1]);
            symmoperations.append(np.copy(symatrix[t,:,:]));
        if len(listline)==9:
          for t in range(1):
            for m in range(3):
              for n in range(3):
                symatrix[t][m][n]=float(listline[t*3*3+3*m+n]);
            symmoperations.append(np.copy(symatrix[t,:,:]));
        k=k+1;
      firstfind=0;
  return symmoperations;
if __name__ == '__main__':
  natom=12;
  masslist=[];
  for i in range(2):
    masslist.append(178.49);
  for i in range(2,4):
    masslist.append(91.224);
  for i in range(4,12):
    masslist.append(15.9999);
  dyn=obtaindyn(natom,"dfpt0.abo");
  for i in range(natom):
    for j in range(natom):
      for m in range(3):
        for n in range(3):
          dyn[i][j][m][n]=dyn[i][j][m][n]/np.sqrt(masslist[i]*masslist[j]);
  dyn2D=np.zeros((3*natom,3*natom));
  for i in range(natom):
    for j in range(natom):
      dyn2D[(i*3):((i+1)*3),(j*3):((j+1)*3)]=np.copy(dyn[i,j,0:3,0:3]);
  w,v=np.linalg.eig(dyn2D);
  axis=readaxis("scfpositive");
  Hposition=readposition(axis,"dfpt0",natom);
  Pposition=readposition(axis,"scfpositive",natom);
  Tposition=readposition(axis,"scfTphase",natom);
  washback(Pposition,natom);
  washback(Tposition,natom);
  washback(Hposition,natom);
  dpT=diffposition(Tposition,Hposition,axis,natom);
  dpP=diffposition(Pposition,Hposition,axis,natom);
  dp=dpP-dpT;
  pchange=dp.reshape((1,3*natom));
  coeff=np.matmul(pchange,np.linalg.inv(v));
  re=np.zeros(3*natom);
else:
  print("File Imported")
  '''
  for i in range(length):
    temp="";
    for m in range(3):
      for n in range(3):
        temp=temp+" "+str(symop[i][m][n]);
    print(temp)
  '''
