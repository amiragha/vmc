import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft
from numpy import pi
from itertools import product

import h5py

#ilename = "vmc3046,kagomestripLC,Lx=32open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"
filename = "vmc9295,kagomestripLC,Lx=40open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"

h5file = h5py.File(filename, 'r')
corr = h5file["correlations"][:]
corr *= 3.
print(corr)
h5file.close()

def index_of(m,n,N):
    if m > n :
        return index_of(n,m,N)
    else:
        return ((2*N-m)*(m-1))//2 + (n-m) - 1

def extract_correlations(corr, N, pattern = "MM"):
    if pattern == "MM":
        correlations = np.eye(N) * 0.75
        for i,j in product(range(N), range(N)) :
            if i==j :
                continue

            i_middle = 3*i + 2
            j_middle = 3*j + 2
            correlations[i-1, j-1] = corr[index_of(i_middle, j_middle, 3*N)]

        ys = np.zeros(N//2+1)
        avg_sfactor = 0
        for i in range(N):
            ys = rfft(np.roll(correlations[:,i],-i))
            avg_sfactor += np.real(ys)

        xs = np.linspace(0, pi, N//2+1)
        return xs, avg_sfactor/N

    if pattern == "lowerleg":
        bonds = np.zeros(N)
        for i in range(N-1):
            i_low = 3*i + 1
            j_low = 3*(i+1) + 1
            bonds[i+1] = corr[index_of(i_low, j_low, 3*N)]

        ys = np.real(rfft(bonds))
        xs = np.linspace(0, pi, N//2+1)

        return xs[1:], ys[1:]
N = 40

fig = plt.figure()
ax = fig.add_subplot(111)

#xs, ys = extract_correlations(corr, N, "MM")
xs, ys = extract_correlations(corr, N, "lowerleg")
ax.plot(xs, ys, "b.")
ax.plot(xs, ys, "b:")

ax.axvline(pi/2, c='k', ls=':')
ax.set_xlim(0,pi)
ax.set_xticks([0, pi])
ax.set_xticklabels(["$0$",r"$\pi$"])

fig.savefig("python_plot.pdf")
