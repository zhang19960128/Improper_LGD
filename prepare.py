import numpy as np
from itertools import combinations
modeid=[11,13,15,19,21,22,29];
comb=combinations(modeid,1);
for i in comb:
  print(list(i))
