import numpy as np
import matplotlib.pyplot as plt
from banded import banded
import scipy.linalg as linalg

#Schrodinger's Eq.: -h^2/2M * ∂^2ψ/∂x^2 = ihbar*∂ψ/∂t
#two walls: at x=0,L

def get_A(L, m, hbar, h, N, a):
    a1 = 1 + (h*1j*hbar/(2*m*a**2))
    a2 = -h*1j*hbar/(4*m*a**2)
    #define A
    A = np.empty((3, N), complex)
    A[0, :] = a2
    A[1, :] = a1
    A[2, :] = a2
    print(A)
    return A

def psi0(x, L):
    '''initial wave function'''
    sigma = 1e-10 #m (a measure of uncertainty in position)
    kappa = 5e10 #1/m (like a wavenumber)
    x0 = L/2
    return np.exp(1j*kappa*x) * np.exp(-(x-x0)**2/(2*sigma**2))

def run(T):
    #define quantities:
    L = 1e-8 #m
    m = 9.109e-31 #kg
    hbar = 1.0546e-36 #Js
    h = 1e-17 #time step #algorithm is numerically unstable for "large" timesteps
    N = 1000 #total number of steps along x
    a = L / N #spacing

    #initialize
    xlist = np.linspace(0, L, N+1)
    psi = psi0(xlist, L)
    #apply BCs
    psi[[0,0]], psi[[0,N]] = 0, 0

    b1 = 1 - (h*1j*hbar/(2*m*a**2))
    b2 = h*1j*hbar/(4*m*a**2)
    A = get_A(L, m, hbar, h, N, a)

    psi_dict = {}

    for t in range(T):
        psi_dict[t] = psi
        v = b1*psi[1:N] + b2*(psi[2:N+1] + psi[:N-1]) #from p441
        psi[1:N] = banded(A, v, 1, 1) #function to solve banded linear system of eq.
            #can also use scipy.linalg.solve_banded
        psi_dict[t] = psi
        print(np.mean(np.real(banded(A, v, 1, 1))))
    return psi_dict

def plotter(psi_dict, index):
    fig, ax = plt.subplots()
    psi = psi_dict[index]
    h = 1e-17 #time step
    L = 1e-8 #m
    N = 1000 #total number of steps
    xlist = np.linspace(0, L, N+1)
    mag = np.sqrt(np.imag(psi)**2 + np.real(psi)**2)
    re_wave = np.real(psi)
    t = h*index
    ax.plot(xlist, re_wave, color = 'navy', alpha = .5, label=f"t = {t}s")
    ax.axvline(x=0, ymin=-1, ymax=1, color = 'k', linestyle = "--")
    ax.axvline(x=L, ymin=-1, ymax=1, color = 'k', linestyle = "--")
    ax.set_xlabel(r'$x$ [m]')
    ax.set_ylabel(r'Re$\{\Psi\}$')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f're_wavefunction_at_time={t}.png', dpi=200)
    #the wave should spread, propagate to the right,
    #and reflect off of the RHS wall.

if __name__ == "__main__":
    T = 50000
    psi_dict = run(T)
    plotter(psi_dict, T-1)