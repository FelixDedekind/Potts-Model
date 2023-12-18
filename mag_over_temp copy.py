import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

df = pd.read_csv("mag_over_temp.txt",delim_whitespace=True,header=None)
T = df.to_numpy()[:,0]
mag = abs(df.to_numpy()[:,1])

plt.scatter(T,mag)
plt.grid()
plt.xlabel("temperature (kB = 1)")
plt.ylabel("magnetization")
plt.xlim(T[0],T[-1])
plt.tight_layout()
plt.show()