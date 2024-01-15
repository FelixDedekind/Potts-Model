import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("lcs_over_temp.txt",delim_whitespace=True,header=None)
T = df.to_numpy()[:,0]
lcs = df.to_numpy()[:,1]
deltalcs = df.to_numpy()[:,2]

plt.errorbar(T,lcs,yerr=deltalcs,fmt='.')
plt.grid()
plt.xlabel("temperature (kB = 1)")
plt.ylabel("largest cluster size")
plt.xlim(T[0],T[-1])
plt.tight_layout()
plt.show()