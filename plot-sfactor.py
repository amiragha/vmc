import numpy as np
import matplotlib.pyplot as plt
from numpy.fft import rfft
from numpy import pi
from itertools import product

import h5py

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

        np.savetxt("vmc,kagomestripLC,Lx=64open,mu={:.1f},tcr={:.1f}.dat".format(mu, tc), ys)
        return xs[1:], ys[1:]

#filename = "vmc1526,kagomestripLC,Lx=32open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"
#filename = "vmc7921,kagomestripLC,Lx=32open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"
#filename = "vmc9295,kagomestripLC,Lx=40open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"

# tclist = [0.6, 1.0, 1.4, 1.8]
# mulist = np.arange(0.1, 5.05, 0.1)

tcmulist = [
    [0.5, -1.2],
    [0.5, -1.4],
    [0.5, -1.7],
    [0.5, -2.1],
    [0.6, -1.3],
    [0.6, -1.5],
    [0.6, -1.9],
    [0.6, -2.5],
    [0.7, -1.4],
    [0.7, -1.7],
    [0.7, -2.1],
    [0.7, -3.1],
    [0.8, -1.5],
    [0.8, -2.0],
    [0.8, -2.3],
    [0.8, -3.6],
    [0.9, -1.6],
    [0.9, -2.1],
    [0.9, -2.7],
    [0.9, -4.1],
    [1.0, -1.8],
    [1.0, -2.4],
    [1.0, -3.1],
    [1.0, -4.8],
    [1.1, -1.9],
    [1.1, -2.6],
    [1.1, -3.4],
    [1.1, -5.6],
    [1.2, -2.1],
    [1.2, -2.9],
    [1.2, -3.9],
    [1.2, -6.4],
    [1.3, -2.3],
    [1.3, -3.2],
    [1.3, -4.3],
    [1.3, -7.0]
]

#for tc, mu in product(tclist, mulist):
for tc, mu in tcmulist:
    filename = "vmc,kagomestripLC,Lx=64open,mu={:.1f},tcr={:.1f},tcphi=0.0.h5".format(mu, tc)
    h5file = h5py.File(filename, 'r')
    corr = h5file["correlations"][:]
    corr *= 3.
    #print(corr)
    h5file.close()

    N = 64

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

    fig.savefig("mu={:.1f},tc={:.1f}.png".format(mu, tc))
    plt.close(fig)
