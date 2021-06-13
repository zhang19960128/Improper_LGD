import numpy as np
import os
def find_nth(string,substring,n):
  if(np.abs(n-1)<1e-6):
    return string.find(substring)
  else:
    previous=find_nth(string,substring,n-1);
    if np.abs(previous+1) < 1e-6:
      return -1;
    else:
      return string.find(substring,previous+1);
def find_symop(line):
  oplist=[];
  for i in range(1,10):
    if find_nth(line,"(",i)!=-1:
      indexone=find_nth(line,"(",i);
      indextwo=find_nth(line,")",i);
      substring=line[(indexone+1):indextwo].replace("|"," ").replace(','," ");
      oplist.append(substring);
  return oplist
def allop(lines):
  optotal=[];
  for i in range(len(lines)):
    oplist=find_symop(lines[i]);
    optotal=optotal+oplist;
  return optotal
def ired(irename):
  iretab=open(irename,'r');
  lines=iretab.readlines();
  for i in range(len(lines)):
    if lines[i].find("****Irrep")!=-1:
      start=i;
  repname=[];
  for i in range(start+1,len(lines)-1):
    repname.append(lines[i].replace("\n",""));
  return repname;
def elements(filename):
  ftab=open(filename,'r');
  lines=ftab.readlines();
  opall=allop(lines[6:]);
  return opall;
def writecombination(irep,op,filename):
  f=open(filename,'w');
  f.write("VALUE PARENT 225\n")
  f.write("VALUE IRREP "+irep+"\n")
  f.write("SHOW MATRIX\n")
  f.write("VALUE ELEMENTS "+op+"\n")
  f.write("DISPLAY IRREP\n")
  f.write("QUIT\n");
  f.close();
def extractmatrix(filename):
  f=open(filename,'r');
  lines=f.readlines();
  for i in range(len(lines)):
    if lines[i].find("Element")!=-1 and lines[i].find("Matrix")!=-1:
      startindex=i;
  dim=len(lines)-1-startindex-2;
  opmatrix=np.zeros((dim,dim));
  for i in range(dim):
    for j in range(dim):
      opmatrix[i][j-dim]=float(lines[startindex+1+i].split()[j-dim]);
  return opmatrix;
def ireduciblerepre():
  folder="IREELE";
  totalelements=elements(str(folder)+"/elements.out");
  ireduce=ired(str(folder)+"/ire.out");
  iredeleop=[];
  for i in range(len(ireduce)):
    ire=[];
    for j in range(len(totalelements)):
      writecombination(ireduce[i],totalelements[j],folder+"/"+"IRE"+str(i)+"ELE"+str(j));
  #    cmd="./iso < IRE"+str(i)+"ELE"+str(j)+" >"+"./IREELE/IRE"+str(i)+"ELE"+str(j)+"out";
  #    os.system(cmd);
      opmatrix=extractmatrix("./"+str(folder)+"/IRE"+str(i)+"ELE"+str(j)+"out");
      ire.append(opmatrix);
    iredeleop.append(ire);
  return [iredeleop,ireduce,totalelements]
def matrixrepresentation():
  f=open("./ISODATA/data_space.txt",'r');
  lines=f.readlines();
  pgroupstart=1;
  for i in range(len(lines)):
    if lines[i].find("point_op_label_stokes")!=-1:
      pgroupend=i;
    if lines[i].find("ipoint_op")!=-1:
      matrixgroupstart=59+1;
    if lines[i].find("ipoint_op_mlt")!=-1:
      matrixgroupend=i;
  labelop=[];
  for i in range(pgroupstart,pgroupend):
    labelop=labelop+lines[i].replace('"',' ').split();
  matrix=[];
  print(matrixgroupstart,matrixgroupend)
  for i in range(matrixgroupstart,matrixgroupend):
    matrix=matrix+lines[i].split();
    print(len(lines[i].split()))
  for i in range(len(matrix)):
    matrix[i]=float(matrix[i]);
  print('length of symop is',len(labelop),"length of matrix is: ",len(matrix)/9)
  matrixlist=[];
  dictionary=dict();
  for i in range(len(labelop)):
    mat=np.zeros((3,3));
    for m in range(3):
      for n in range(3):
        mat[m][n]=matrix[i*3*3+m*3+n];
    dictionary[labelop[i]]=mat;
  return dictionary;
[iredeleop,ireduce,elements]=ireduciblerepre();
dictionary=matrixrepresentation();
for i in range(len(dictionary)):
  print(dictionary[i])
