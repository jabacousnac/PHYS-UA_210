import numpy as np
import matplotlib.pyplot as plt

def radioactive():
    #Tl208 -> Pb208, with half life of 3.053 min

    N0 = 1000
    lambda_Tl208 = 3.053*60 #in seconds
    mu = np.log(2)/lambda_Tl208 #from p459 in the book
    
    z = np.random.random(N0) #generate 1000 random numbers
    t = -1/mu * np.log(1-z)
    t = np.sort(t)
    survived = N0-np.arange(1, N0+1)
    
    plt.figure()
    plt.plot(t/60, survived)
    plt.ylabel('N')
    plt.xlabel('t [min]')
    plt.savefig('radioactive02.png')

if __name__ == "__main__":
    radioactive()
