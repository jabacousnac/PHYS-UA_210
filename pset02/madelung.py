import numpy as np
import matplotlib.pyplot as plt
import itertools
import timeit

#the madelung constant (Ex: 2.9)

def with_for_loop(L):
    start = timeit.default_timer()
    M = 0
    L = int(L)
    the_range = range(-L, L)
    for i in the_range:
        for j in the_range:
            for k in the_range:
                if (i == 0 and j == 0 and k == 0):
                    continue
                else:
                    sgn = (-1)**(abs(i+j+k)%2)
                    M -= sgn/np.sqrt(i**2 + j**2 + k**2)
    dt = timeit.default_timer() - start
    print('L = {}\t'.format(L) + 'M = {}'.format(M))
    print(str(dt) + 's\n')

def without_for_loop(L):
    start = timeit.default_timer()
    L = int(L)
    the_range = range(-L, L)
    i, j, k = np.meshgrid(the_range, the_range, the_range)
    mask = ~(i==0) | ~(j==0) | ~(k==0)
    print(mask)
    sgn = (-1)**(np.abs(i+j+k)%2)
    M = -np.sum(sgn[mask] / np.sqrt(i[mask]**2 + j[mask]**2 + k[mask]**2))
    dt = timeit.default_timer() - start
    print('L = {}\t'.format(L) + 'M = {}'.format(M))
    print(str(dt) + 's\n')

if __name__ == "__main__":
    with_for_loop(1e1)
    without_for_loop(1e1)
