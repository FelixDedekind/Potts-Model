import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("energy_over_temp.txt",delim_whitespace=True,header=None)
T = df.to_numpy()[:,0]
e = df.to_numpy()[:,1]

plt.scatter(T,e)
plt.grid()
plt.xlabel("temperature (kB = 1)")
plt.ylabel("energy")
plt.xlim(T[0],T[-1])
plt.tight_layout()
plt.show()