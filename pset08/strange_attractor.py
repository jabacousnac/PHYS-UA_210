import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

#Lorenz Eq.
#dx/dt = sigma(y-x), dy/dt = rx-y-xz, dz/dt = xy-bz; sigma, r, b are constants.
def lorenz(t, X, sigma, r, b):
    x, y, z = X
    dxdt = sigma*(y-x)
    dydt = r*x - y - x*z
    dzdt = x*y - b*z
    return dxdt, dydt, dzdt #this is the integrand

def solve_lorenz(sigma, r, b, v):
    #set up IVP.
    tmax = 50
    x0, y0, z0 = 0, 1, 0
    ans = solve_ivp(lorenz, (0, tmax), (x0, y0, z0), args = (sigma, r, b), dense_output = True) #first arg is the function, followed by time interval, then initial values, then values for arguments
    #plot y v. t and z v. x:
    t = np.linspace(0, tmax, 10000)
    x, y, z = ans.sol(t)
    fig, axs = plt.subplots(2,1)
    #y v. t
    ax = axs[0]
    ax.set_xlabel('t')
    ax.set_ylabel('y')
    ax.plot(t, y, 'r-', label = r'$\sigma = {},\,$'.format(sigma) + '$r = {},\,$'.format(r) + '$b = {}$'.format(np.around(b,3)))
    #z v. x
    ax = axs[1]
    ax.set_xlabel('x')
    ax.set_ylabel('z')
    p = ax.scatter(x, z, marker = 'o', s = 4, c=t, alpha = .4, cmap = 'magma', label = r'$\sigma = {},\,$'.format(sigma) + '$r = {},\,$'.format(r) + '$b = {}$'.format(np.around(b,3)))
    plt.colorbar(p)
    for ax in axs:
        ax.legend()
    plt.tight_layout()
    plt.savefig(f'strange_attractor{str(v).zfill(2)}.png', dpi=200)

if __name__ == "__main__":
    solve_lorenz(10, 28, 8/3, 1)
