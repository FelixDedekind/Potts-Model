import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("mag_over_temp.txt",delim_whitespace=True,header=None)
T = df.to_numpy()[:,0]
mag = df.to_numpy()[:,1]
deltamag = df.to_numpy()[:,2]

plt.errorbar(T,mag,yerr=deltamag,fmt='.')
plt.grid()
plt.xlabel("temperature (kB = 1)")
plt.ylabel("magnetization")
plt.xlim(T[0],T[-1])
plt.tight_layout()
plt.show()