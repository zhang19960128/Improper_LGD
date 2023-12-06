import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt("Enthalpy.dat")
plt.plot(np.arange(len(data)),(data-np.min(data))*13.59*2)
plt.xlabel("PATH")
plt.ylabel("Enthalpy(eV)")
plt.show()
