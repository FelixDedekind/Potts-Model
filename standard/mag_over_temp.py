import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib

font = {'size'   : 18}
matplotlib.rc('font', **font)

df = pd.read_csv("q=4_closeup/mag_over_temp.txt",delim_whitespace=True,header=None)
T = df.to_numpy()[:,0]
mag = df.to_numpy()[:,1]
deltamag = df.to_numpy()[:,2]

#plt.errorbar(T,mag,yerr=deltamag,fmt='.')
plt.scatter(T,mag)
plt.grid()
plt.xlabel(r"$\mathrm{k_BT}/\mathrm{\epsilon}$")
plt.ylabel("normalized magnetization")
plt.xlim(T[0],T[-1])
plt.ylim(0,1)
plt.tight_layout()
plt.show()