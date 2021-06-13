import numpy as np
from tabulate import tabulate
def printoscreencharater(operations,outmatrix):
  sp=np.shape(outmatrix);
  total=[];
  head=[];
  for i in range(sp[0]):
    line=[];
    line.append(operations[i]);
    for j in range(sp[1]):
      line.append("{0:3.2f}".format(outmatrix[i][j]));
    total.append(line);
  head.append("Symm Ops");
  for i in range(sp[1]):
    head.append("MO"+str(i));
  print (tabulate(total, headers=head))
