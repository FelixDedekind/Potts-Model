import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib

font = {'size'   : 18}
matplotlib.rc('font', **font)

df = pd.read_csv("q=4_closeup/lcs_over_temp.txt",delim_whitespace=True,header=None)
T = df.to_numpy()[:,0]
lcs = df.to_numpy()[:,1]
deltalcs = df.to_numpy()[:,2]

#plt.errorbar(T,lcs,yerr=deltalcs,fmt='.')
plt.scatter(T,lcs)
plt.grid()
plt.xlabel(r"$\mathrm{k_BT}/\mathrm{\epsilon}$")
plt.ylabel("normalized largest cluster size")
plt.xlim(T[0]-0.02,T[-1]+0.02)
plt.tight_layout()
plt.show()