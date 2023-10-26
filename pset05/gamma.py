import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quadrature

def integrand(a, x):
    return x**(a-1)*np.exp(-x)

def plotter():
    fig, ax = plt.subplots()
    xlist = np.linspace(0, 5, 100)
    y2 = [integrand(2, x) for x in xlist]
    y3 = [integrand(3, x) for x in xlist]
    y4 = [integrand(4, x) for x in xlist]
    ax.plot(xlist, y2, 'r-')
    ax.plot(xlist, y3, 'g-')
    ax.plot(xlist, y4, 'b-')
    ax.set_xlabel('x')
    ax.set_ylabel('integrand')
    plt.savefig('integrand.png', dpi=96)
    
    
def integrate(a):
    f = lambda z: np.exp((a-1)*np.log(z*(a-1)/(1-z)) - z*(a-1)/(1-z)) * (a-1)/(1-z)**2
    res = quadrature(f, 0., 1.)
    print(res)
    
    
if __name__ == "__main__":
    plotter()
    print(integrate(3/2)[0])
