import numpy as np
import matplotlib.pyplot as plt
import timeit

def with_for_loops(N):
    start = timeit.default_timer()
    A, B, C = np.zeros([N, N], float), np.zeros([N, N], float), np.zeros([N, N], float)
    for i in range(N):
        for j in range(N):
            for k in range(N):
                C[i,j] += A[i,k]*B[k,j]
    dt = timeit.default_timer() - start
    return C, dt

def without_for_loops(N):
    start = timeit.default_timer()
    A, B = np.zeros([N, N], float), np.zeros([N, N], float)
    dt = timeit.default_timer() - start
    return np.dot(A, B), dt

if __name__ == "__main__":
    plt.figure()
    Nlist = [10, 30, 75, 150, 300, 500]
    dt = []
    for N in Nlist:
        print('we are at N = {}'.format(N))
        dt.append(with_for_loops(N)[1]) #change method here
        print(N, dt)
    plt.scatter(np.log10(Nlist), np.log10(dt))
    m, b = np.polyfit(np.log10(Nlist), np.log10(dt), 1)
    xfit = np.linspace(0, 3, 100)
    yfit = xfit*m + b
    plt.plot(xfit, yfit, 'k-', label = r'slope = {}'.format(np.around(m,3)))
    plt.xlabel('$log_{10}\,(N)$')
    plt.ylabel('$log_{10}\,(dt [s])$')
    plt.legend()
    plt.tight_layout()
    plt.savefig('matrix_mult_with_for_loops.png') #change filename, correspondingly
