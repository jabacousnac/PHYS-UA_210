import numpy as np
import matplotlib.pyplot as plt

def mandelbrot(c, N):
    dim = c.shape
    grid = np.resize(np.array(0,), dim)
    z = np.zeros(dim, dtype=np.complex64)
    for i in range(N):
        z = z**2 + c
        done = np.greater(abs(z), 2.0)
        c = np.where(done, 0+0j, c)
        z = np.where(done, 0+0j, z)
        grid = np.where(done, i, grid)
    return grid

def draw(N):
    x = np.linspace(-2, 2, 100)
    y = np.linspace(-2, 2, 100)
    c = x[:, np.newaxis] + y*1j
    grid = mandelbrot(c,N)
    
    #plot:
    fig, ax = plt.subplots()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.imshow(grid.T,cmap='jet')
    plt.savefig('mandelbrot2.png')
    return

    
if __name__ == "__main__":
    draw(100)
