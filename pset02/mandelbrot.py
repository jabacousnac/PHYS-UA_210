import numpy as np
import matplotlib.pyplot as plt

def mandelbrot(N):
    x = np.linspace(-2, 2, N)
    y = np.linspace(-2, 2, N)
    c = x[:, np.newaxis] + 1j*y
    #we have generated c above. Now we need to iterate
    grid = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            z = 0
            this_c = c[i, j]
            n = 0
            while n < 100:
                n += 1
                z = np.real(z**2) + np.real(this_c) + 1j*(np.imag(z**2) + np.imag(this_c))
                if check_iteration(z):
                    grid[j, i] = n
                    n = 101
    print(grid)
    #plotting
    plt.figure()
    plt.xlabel('x')
    plt.ylabel('y')
    p = plt.imshow(grid, cmap = 'jet')
    q = plt.colorbar(p)
    q.set_label('n')
    plt.savefig('mandelbrot_n={}.png'.format(N))
    
def check_iteration(z):
    if abs(z) > 2:
        return True
    else:
        return False

if __name__ == "__main__":
    mandelbrot(1000)
