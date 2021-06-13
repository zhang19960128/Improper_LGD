import numpy as np
f=open("POSCAR",'r');
lines=f.readlines();
outlist=[];
for i in range(len(lines)-1,7,-1):
  line=lines[i].split();
  temp="";
  for j in range(3):
    re=float(line[j])-float(lines[8].split()[j]);
    if re < 0:
      re= re +1;
    elif re > 1:
      re= re -1;
    else:
      pass;
    temp=temp+"{0:10.7f}".format(re)+" ";
  temp=temp+line[3];
  outlist.append(temp);
outlist.reverse();
for i in range(0,8):
  print(lines[i].rstrip())
for i in range(len(outlist)):
  print(outlist[i])
