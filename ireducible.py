import numpy as np
import re
def readirreducible():
  webfile="website.txt";
  f=open(webfile,'r');
  lines=f.readlines();
  nameofirep=[];
  for i in range(len(lines)):
    if lines[i].find("Number of elements")!=-1:
      numofsym=int(lines[i].split()[-1]);
    if lines[i].find("allowed irreps")!=-1:
      numofrep=int(lines[i].split()[-3]);
  print(numofsym,numofrep);
  dimension=[];
  opmatrixinrep=[];
  for i in range(len(lines)):
    if lines[i].find("Irrep")!=-1 and lines[i].find("dimension")!=-1:
      nameofirep.append(lines[i].split()[1]);
      opcount=0;
      dim=int(lines[i].split()[-1]);
      dimension.append(dim);
      opmatrix=[];
      linenum=i+2;
      while opcount!=numofsym:
        line=lines[linenum];
        line=line.replace("(","");
        line=line.replace(")","");
        line=line.replace(","," ");
        symop=len(line.split())/2/dim;
        symop=int(symop);
        opcount=opcount+symop;
        linecollect=[];
        for j in range(dim):
          line=lines[linenum+j];
          line=line.replace("(","");
          line=line.replace(")","");
          line=line.replace(","," ");
          line=line.replace("\n"," ");
          linecollect.append([float(i) for i in line.split()]);
        for j in range(symop):
          work=np.zeros((dim,2*dim));
          symop=np.zeros((dim,dim),dtype=complex);
          for m in range(dim):
            for n in range(2*dim):
              work[m][n]=linecollect[m][j*2*dim+n];
          for m in range(dim):
            for n in range(dim):
              symop[m][n]=work[m][2*n]*(np.cos(work[m][2*n+1]/180*3.141592653)+complex(0,1)*np.sin(work[m][2*n+1]/180*3.141592653));
          opmatrix.append(symop);
        linenum=linenum+dim+2;
      opmatrixinrep.append(opmatrix);
  opmatrix=[];
  for i in range(len(lines)):
    if lines[i].find("elements as translation coset representatives :")!=-1:
      startindex=i+3;
      symcount=0;
      while symcount != numofsym:
        numop=len(lines[startindex].split())/4;
        print(numop)
        for j in range(int(numop)):
          symop=np.zeros((3,3));
          for m in range(3):
            for n in range(3):
              symop[m][n]=float(lines[m+startindex].split()[j*4+n]);
          opmatrix.append(symop);
          symcount=symcount+1;
        startindex=startindex+5;
  return([opmatrixinrep,opmatrix,nameofirep])
