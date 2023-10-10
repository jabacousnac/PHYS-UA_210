import numpy as np
import matplotlib.pyplot as plt
from gaussxw import gaussxw

def func(a, x):
    return 1/np.sqrt(a**4 - x**4)

def integrate(a):
    N = 20
    xmin = 0
    xmax = a
    x, w = gaussxw(N)
    xp = .5*(xmax-xmin)*x + .5*(xmin+xmax)
    wp = .5*(xmax-xmin)*w
    the_sum = 0.
    for i in range(N):
        the_sum += wp[i]*func(a, xp[i])
    return np.sqrt(8)*the_sum

def plotter():
    y_list = []
    a_list = np.linspace(0, 2, 50)[1:]
    for a in a_list:
        y_list.append(integrate(a))
    fig, ax = plt.subplots()
    ax.grid(alpha = .5)
    ax.set_xlabel(r'$a$')
    ax.set_ylabel(r'$T$')
    ax.scatter(a_list, y_list, marker = 'o', color='hotpink')
    plt.savefig('anharmonic.png')
        
if __name__ == "__main__":
    plotter()
    
    
    
