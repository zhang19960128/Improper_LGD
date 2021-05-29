import analyzePH
import numpy as np
import math
import cmath
import sys
import check
import IOmode
import inputfiles
import primitive
import sequence
import obtainmode
from mpi4py import MPI
import sciconst
import groupmatch
import extend
#Note that symmetry operation also move the atoms#
comm=MPI.COMM_WORLD
size=comm.Get_size();
rank=comm.Get_rank();
[primitiveM,wM,vM]=primitive.obtainprimitivemode(inputfiles.natomPM,inputfiles.dfptinPM,inputfiles.dfptoutPM,inputfiles.modePMname);
[primitiveGa,wGa,vGa]=primitive.obtainprimitivemode(inputfiles.natomPGa,inputfiles.dfptinGa,inputfiles.dfptoutGa,inputfiles.modePGaname);
[masslist,namelist]=sequence.sequence(inputfiles.natomExp);
axis=analyzePH.readaxis(inputfiles.dfptinExp);
[Dw,Dv]=obtainmode.obtainmode(inputfiles.natomExp,masslist,inputfiles.dfptinExp,inputfiles.dfptoutExp);
ExpandPosition=analyzePH.readposition(axis,inputfiles.dfptinExp,inputfiles.natomExp);
EvM=extend.modemap(ExpandPosition,masslist,axis,primitiveM,vM.real,-1);
EvGamma=extend.modemap(ExpandPosition,masslist,axis,primitiveGa,vGa.real,1)
if rank==0:
  IOmode.printmode(ExpandPosition,namelist,axis,EvM,inputfiles.modeExpMname);
  IOmode.printmode(ExpandPosition,namelist,axis,EvGamma,inputfiles.modeExpGaname);
  for i in range(len(wGa)):
    IOmode.printmode(ExpandPosition,namelist,axis,[EvGamma[i]],inputfiles.modeExpGaname+"{0:6.2f}".format(cmath.sqrt(wGa[i])*np.sqrt(sciconst.Ha/sciconst.bohr/sciconst.bohr/sciconst.aumass)*10**-12*33.35641/2/3.141592653));
symop=analyzePH.readsymmetry(inputfiles.dfptoutExp);
length=len(symop)
localsymop=symop[rank:length:size];
localTlist=check.check(localsymop,axis,masslist,ExpandPosition);
print(localsymop[3])
[matchcoeff,matchmode]=groupmatch.groupmatchoperation(wGa,localsymop[3],localTlist[3],EvGamma,axis);
print(len(matchmode),len(EvGamma))
IOmode.printmode(ExpandPosition,namelist,axis,matchmode,"./MATCH/MATCHMODE");
