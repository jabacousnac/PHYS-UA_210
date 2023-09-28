import numpy as np
import matplotlib.pyplot as plt
from random import random

def run():

    N0 = int(1e4)
    dt = 1
    tmax = 5400 #90 min
    tpoints = np.arange(0., tmax, dt)
    
    Bi213, Tl209, Pb209, Bi209 = N0, 0, 0, 0
    Bi213_points, Tl209_points, Pb209_points, Bi209_points = [N0], [0], [0], [0]
    lambda_Bi213, lambda_Tl209, lambda_Pb209 = 46*60, 2.2*60, 3.3*60
    
    #start the main loop
    for t in tpoints:

        if t%1000 == 0:
            print("we are at t = {}s".format(t))
            
        #do the Pb209 -> Bi209 decay (half life is 3.3 min)
        p = 1 - 2**(-dt/lambda_Pb209)
        decay = 0
        for i in range(Pb209):
            if random() < p:
                decay += 1
        Pb209 -= decay
        Bi209 += decay

        #do the Tl209 -> Pb209 decay (half life is 2.2 min)
        p = 1 - 2**(-dt/lambda_Tl209)
        decay = 0
        for i in range(Tl209):
            if random() < p:
                decay += 1
        Tl209 -= decay
        Pb209 += decay

        #finally, do the Bi213 -> Pb209 decay (half life is 46 min), but also the\
            #Bi213 -> Tl209 pathway
        p = (1 - 2**(-dt/lambda_Bi213))
        decay = 0
        for i in range(Bi213):
            q = random()
            if q > p: #it does not decay
                continue
            else: #it decays, and we have to choose a pathway
                decay += 1
                if q/p < .9791:
                    Pb209 += 1
                else:
                    Tl209 += 1
        Bi213 -= decay            
            
        #append the lists 
        Bi213_points.append(Bi213)
        Tl209_points.append(Tl209)
        Pb209_points.append(Pb209)
        Bi209_points.append(Bi209)

    plot(Bi213_points, Tl209_points, Pb209_points, Bi209_points)
    
    return


def plot(Y0, Y1, Y2, Y3):

    plt.figure()
    t = np.arange(len(Y0))/60
    plt.plot(t, Y0, 'k-', label = 'Bi-213')
    plt.plot(t, Y1, 'r-', label = 'Tl-209')
    plt.plot(t, Y2, 'g-', label = 'Pb-209')
    plt.plot(t, Y3, 'b-', label = 'Bi-209')
    
    plt.xlabel('t [min]')
    plt.ylabel('N')
    plt.legend()

    plt.tight_layout()
    plt.savefig('radioactive.png')

if __name__ == "__main__":
    run()
