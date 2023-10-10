#the integrand: x^4 - 2x + 1, limits: x=0...2
#exact answer is 4.4

import numpy as np

def trapz(a, b, n):
    h = (b-a)/n
    xarray = np.arange(a+h, b, h)
    print(h, xarray)
    sums = np.array([f(x) for x in xarray])
    integral = h * (.5*(f(a) + f(b)) + np.sum(sums))
    return integral
    
def f(x):
    return x**4 - 2*x + 1

if __name__ == "__main__":
    I1 = trapz(0, 2, 10)
    I2 = trapz(0, 2, 20)
    print('n = {}'.format(10) + '\t res = {}'.format(I1))
    print('n = {}'.format(20) + '\t res = {}'.format(I2))
    err = 1/3 * (I2 - I1)
    print('err = {}'.format(abs(err)))
