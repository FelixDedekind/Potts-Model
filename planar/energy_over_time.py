import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("energy_over_time.txt",delim_whitespace=True,header=None)
t = df.to_numpy()[:,0]
e = df.to_numpy()[:,1]

plt.plot(t,e)
plt.grid()
plt.xlabel("mc timesteps")
plt.ylabel("energy")
plt.xlim(t[0],t[-1])
plt.tight_layout()
plt.show()