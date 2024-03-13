import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib

font = {'size'   : 18}
matplotlib.rc('font', **font)

q = [2,3,4,5]
Tc = [1.15,1.025,0.925,0.8625]
Tc_error = [0.035,0.02,0.025,0.0125]

qspace = np.linspace(q[0]-1,q[-1]+1,1000)
Tcspace = [1/(np.log(1+np.sqrt(qq))) for qq in qspace]

plt.plot(qspace,Tcspace,label="$\mathrm{log}(1+\sqrt{\mathrm{q}})^{-1}$")
plt.errorbar(q,Tc,yerr=Tc_error,capsize = 5,fmt = '.',label=r"$\mathrm{T_c}$")
plt.xlim(q[0]-0.5,q[-1]+0.5)
plt.grid()
plt.legend()
plt.show()
