import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x*(x-1)

def numerical_fprime(x, delta):
    return (f(x+delta) - f(x))/delta

if __name__ == "__main__":
    print(numerical_fprime(1, 1e-2))
    plt.figure()
    for delta in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14]:
        y = numerical_fprime(1, delta)
        plt.scatter(np.log10(delta), np.log(abs(y-1.)))
    plt.xlabel(r'$log_{10}(\delta)$')
    plt.ylabel(r'$log_{10}$|err| in $df/dx$')
    plt.savefig('calc_derivatives.png')
