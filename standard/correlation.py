import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 
import matplotlib

font = {'size'   : 18}
matplotlib.rc('font', **font)

df = pd.read_csv("correlation.txt",delim_whitespace=True,header=None)
C = df.to_numpy()[:,1]
t = df.to_numpy()[:,0]

plt.scatter(t,C)
plt.grid()
plt.xscale("log")
plt.xlabel("MC timesteps")
plt.ylabel("Correlation")
plt.tight_layout()
plt.show()