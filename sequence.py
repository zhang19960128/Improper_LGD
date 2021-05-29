def sequence(natom):
  masslist=[];
  namelist=[];
  for i in range(int(1/6*natom)):
    masslist.append(178.49);
    namelist.append("Hf")
  for i in range(int(1/6*natom),int(2/6*natom)):
    masslist.append(91.224);
    namelist.append("Zr")
  for i in range(int(2/6*natom),natom):
    masslist.append(15.9999);
    namelist.append("O")
  return([masslist,namelist]);
