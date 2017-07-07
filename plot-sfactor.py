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
        np.savetxt("vmc-ssFT-MM-Lx={}APBC-mu={:.1f},tcr={:.1f}.dat".format(N,mu, tc), avg_sfactor/N)
        return xs, avg_sfactor/N

    elif pattern == "BB":
        correlations = np.eye(N) * 0.75
        for i,j in product(range(N), range(N)) :
            if i==j :
                continue

            i_bottom = 3*i + 1
            j_bottom = 3*j + 1
            correlations[i-1, j-1] = corr[index_of(i_bottom, j_bottom, 3*N)]

        ys = np.zeros(N//2+1)
        avg_sfactor = 0
        for i in range(N):
            ys = rfft(np.roll(correlations[:,i],-i))
            avg_sfactor += np.real(ys)

        xs = np.linspace(0, pi, N//2+1)
        np.savetxt("vmc-ssFT-BB-Lx={}APBC-mu={:.1f},tcr={:.1f}.dat".format(N,mu, tc), avg_sfactor/N)
        return xs, avg_sfactor/N

    elif pattern == "BM":
        correlations = np.zeros((N, N))
        offset_BM = np.array([np.exp(-1.j*2*pi*k*0.5/N) for k in range(0,N//2+1)])
        print (offset_BM)
        for i,j in product(range(N), range(N)) :

            i_bottom = 3*i + 1
            j_middle = 3*j + 2
            correlations[i-1, j-1] = corr[index_of(i_bottom, j_middle, 3*N)]

        ys = np.zeros(N//2+1)
        avg_sfactor = np.zeros(N//2+1)
        for i in range(N):
            ys = rfft(np.roll(correlations[:,i],-i)) * offset_BM
            avg_sfactor += np.real(ys)

        xs = np.linspace(0, pi, N//2+1)
        np.savetxt("vmc-ssFT-BM-Lx={}APBC-mu={:.1f},tcr={:.1f}.dat".format(N,mu, tc), avg_sfactor/N)
        return xs, avg_sfactor/N

    elif pattern == "BT":
        correlations = np.zeros((N, N))
        for i,j in product(range(N), range(N)) :

            i_bottom = 3*i + 1
            j_top = 3*j + 3
            correlations[i-1, j-1] = corr[index_of(i_bottom, j_top, 3*N)]

        ys = np.zeros(N//2+1)
        avg_sfactor = 0
        for i in range(N):
            ys = rfft(np.roll(correlations[:,i],-i))
            avg_sfactor += np.real(ys)

        xs = np.linspace(0, pi, N//2+1)
        np.savetxt("vmc-ssFT-BT-Lx={}APBC-mu={:.1f},tcr={:.1f}.dat".format(N,mu, tc), avg_sfactor/N)
        return xs, avg_sfactor/N

    elif pattern == "lowerleg":
        bonds = np.zeros(N)
        for i in range(N-1):
            i_low = 3*i + 1
            j_low = 3*(i+1) + 1
            bonds[i+1] = corr[index_of(i_low, j_low, 3*N)]

        ys = np.real(rfft(bonds))
        xs = np.linspace(0, pi, N//2+1)

        np.savetxt("vmc-dimerFT-BB-Lx={}APBC-mu={:.1f},tcr={:.1f}.dat".format(N, mu, tc), ys)
        return xs[1:], ys[1:]

#filename = "vmc1526,kagomestripLC,Lx=32open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"
#filename = "vmc7921,kagomestripLC,Lx=32open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"
#filename = "vmc9295,kagomestripLC,Lx=40open,mu=-2.5,tcr=1.0,tcphi=0.0.h5"

tcmulist = [
    [1.0, -2.5]
]

for tc, mu in tcmulist:
    N = 32
    filename = "vmc6722,kagomestripLC,Lx={}periodic,mu={:.1f},tcr={:.1f},tcphi=0.0.h5".format(N, mu, tc)
    h5file = h5py.File(filename, 'r')
    corr = h5file["correlations"][:]
    corr *= 3.
    #print(corr)
    h5file.close()

    fig = plt.figure()
    prows, pcols = 2, 2
    SS = ['BB', 'MM', 'BM', 'BT']
    plot_idx = 0
    for ss in SS:
        plot_idx += 1
        ax = fig.add_subplot(prows, pcols, plot_idx)
        xs, ys = extract_correlations(corr, N, ss)

        ax.plot(xs, ys, "b.")
        ax.plot(xs, ys, "b:")

        ax.axvline(pi/2, c='k', ls=':')
        ax.set_xlim(0,pi)
        ax.set_xticks([0, pi])
        ax.set_xticklabels(["$0$",r"$\pi$"])

    fig.savefig("mu={:.1f},tc={:.1f}.png".format(mu, tc))
    plt.close(fig)
