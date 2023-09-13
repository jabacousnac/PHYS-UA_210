import numpy as np
import matplotlib.pyplot as plt

def gaussian(xlist, mu, sigma):
    A = 1/(np.sqrt(2*np.pi*sigma**2)) #normalization factor
    y = [A*np.exp(-(x-mu)**2/(2*sigma**2)) for x in xlist]
    return y

def plot():
    xlist = np.linspace(-10, 10, 100)
    y = gaussian(xlist, 0, 3)
    #create fig
    plt.figure(0)
    plt.plot(xlist, y, 'ko')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.savefig("gaussian.png")
    plt.show()

if __name__ == "__main__":
    plot()
    
