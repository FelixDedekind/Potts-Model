import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib

font = {'size'   : 18}
matplotlib.rc('font', **font)

q = [2,3,4,5]
Tc = [1.15,1.01,0.925,0.8625]
Tc_error = [0.02,0.02,0.02,0.02]

qspace = np.linspace(q[0]-1,q[-1]+1,1000)
Tcspace = [1/(np.log(1+np.sqrt(qq))) for qq in qspace]

plt.plot(qspace,Tcspace,label="$\mathrm{log}(1+\sqrt{\mathrm{q}})^{-1}$")
plt.errorbar(q,Tc,yerr=Tc_error,capsize = 5,fmt = '.',label=r"simulation")
plt.xlim(q[0]-0.5,q[-1]+0.5)
plt.xlabel("q")
plt.ylabel(r"$\mathrm{k_BT_c/\epsilon}$")
plt.grid()
plt.legend()
plt.show()
