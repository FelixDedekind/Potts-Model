import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

N = [20,25,50, 75,100,125,200,300]
q = [2,3,4,5]
Tcs = np.zeros((len(N),len(q)))

def findmax(array,maxmag):
    tmax = -1
    for cc in range(len(array)):
        if array[len(array)-cc-1] > maxmag/3:
            tmax = cc
            break
    return tmax

for nindex,nn in enumerate(N):
    for qindex,qq in enumerate(q):
        filename = "~/Desktop/temp_potts/Potts/finer/N=" + str(nn) + "/q=" + str(qq) + "/mag_over_temp.txt"
        print(filename)
        df = pd.read_csv(filename, delim_whitespace=True, header=None)
        T = df.to_numpy()[:,0]
        mag = df.to_numpy()[:,1]
        deltamag = df.to_numpy()[:,2]

        maxmag = max(mag)

        deltams = np.array([abs(mag[cc+1] - mag[cc]) for cc in range(len(mag)-1)])
        Tcs[nindex,qindex] = T[findmax(deltams,maxmag)]
        plt.scatter([(T[cc+1]+T[cc])/2 for cc in range(len(T)-1)],deltams)
        plt.plot(T,mag)
        print(T, deltams.argmax())
        plt.show()


