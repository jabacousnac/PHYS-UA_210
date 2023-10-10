import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxw
from scipy.special import roots_hermite

def H(n,x):
    #recursion yay!
    if n==0: #base case 1
        return 1
    elif n==1: #base case 2
        return 2*x
    else:
        return 2*x*(H(n-1, x)) - 2*(n-1)*H(n-2,x)

def factorial(n):
    #can import your function, or you can write your own (like I did)
    #also using recursion (slower than importing from math module)
    if n == 0:
        return 1
    else:
        return n*factorial(n-1)

def psi_func(n, x):
    prefactor = 1/np.sqrt(2**n * factorial(n) * np.sqrt(np.pi))
    return prefactor * np.exp(-(x**2)/2) * H(n,x)
    
def plotting_hermite():
    x_range = np.linspace(-4, 4, 100)
    nlist = [0, 1, 2, 3]
    clist = ['r', 'goldenrod', 'g', 'b']
    #for part b), replace with:
    #nlist, clist, x_range = [30], 'k', np.linspace(-10, 10, 100)
    fig, ax = plt.subplots()
    for i in range(len(nlist)):
        psi = [psi_func(nlist[i], x_range[j]) for j in range(len(x_range))]
        ax.plot(x_range, psi, color = clist[i], label = 'n = {}'.format(nlist[i]))
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$\psi\,(x)$")
    ax.grid(alpha = .5)
    ax.legend()
    plt.tight_layout()
    plt.savefig("wavefunction01.png", dpi=96) #replace it to 'wavefunction02.png' for part b)

def integrand(n, x):
    return abs(psi_func(n,x))**2 * x**2
    
def quant_unc(n, a, b): #let n and the bounds be arguments
    N = 100
    x,w = gaussxw(N)
    xp = .5*(b-a)*x + .5*(b+a)
    wp = .5*(b-a)*w
    the_sum = 0.
    for i in range(N):
        the_sum += wp[i]*integrand(n, xp[i])
    return the_sum

def gauss_hermite(n, a, b):
    N = 10000
    x, w = roots_hermite(N)
    xp = .5*(b-a)*x + .5*(b+a)
    wp = .5*(b-a)*w
    the_sum = 0.
    for i in range(N):
        the_sum += wp[i]*integrand(n, xp[i])
    return the_sum

if __name__ == "__main__":
    #plotting_hermite()
    #print(np.sqrt(quant_unc(5,-10,10))) #should get about 2.345
    print(np.sqrt(gauss_hermite(5, -200, 200)))
