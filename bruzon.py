import numpy as np
def integer(a):
  if a-round(a) < 1e-8:
    return 1;
  else:
    return 0;
def match(qpt1,qpt2):
  matchyes=1;
  for i in range(3):
    matchyes=matchyes*integer(qpt1[i]-qpt2[i]);
  return matchyes;
def groupkpoint(qlist):
  qpointgroup=[];
  for i in range(len(qlist)):
    startmatch=0;
    insertsignal=0;
    while startmatch < len(qpointgroup):
      if match(qlist[i],qpointgroup[startmatch][0]):
        qpointgroup[startmatch].append(qlist[i]);
        insertsignal=1;
        break;
      else:
        insertsignal=0;
      startmatch=startmatch+1;
    if insertsignal==0:
      qpointgroup.append([qlist[i]]);
  reducedkpt=[];
  for i in range(len(qpointgroup)):
    reducedkpt.append(qpointgroup[i][0]);
  return reducedkpt;
basisone=[[0.5,0.5,0],[0,0.5,0.5],[0.5,0,0.5]];
#basisone=[[0.5,0.5,0],[-0.5,0.5,0.0],[0.0,0.0,1.0]];
basistwo=[[1,0,0],[0,1,0],[0,0,1]];
basisone=np.array(basisone);
basistwo=np.array(basistwo);
recipone=np.linalg.inv(basisone).transpose();
reciptwo=np.linalg.inv(basistwo).transpose();
qlist=[];
for m in range(-8,9,1):
  for n in range(-8,9,1):
    for l in range(-8,9,1):
      if m==0:
        scale1=0;
      else:
        scale1=1/m;
      if n==0:
        scale2=0;
      else:
        scale2=1/n;
      if l==0:
        scale3=0;
      else:
        scale3=1/l;
      p1=(scale1*recipone[0]+scale2*recipone[1]+scale3*recipone[2]).dot(reciptwo[0]);
      p2=(scale1*recipone[0]+scale2*recipone[1]+scale3*recipone[2]).dot(reciptwo[1]);
      p3=(scale1*recipone[0]+scale2*recipone[1]+scale3*recipone[2]).dot(reciptwo[2]);
      if np.abs(p1-round(p1)) < 1e-6 and np.abs(p2-round(p2)) < 1e-6 and np.abs(p3-round(p3)) < 1e-6:
        qlist.append([scale1,scale2,scale3]);
print(qlist)
print("------------------------------------")
print(groupkpoint(qlist))
