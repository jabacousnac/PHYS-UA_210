import numpy as np
import matplotlib.pyplot as plt
import scipy as sc

f = lambda x: (x-.3)**2 * np.exp(x) #min at x=.3

def s_quad_interp(a, b, c):
    epsilon = 1e-7 #for numerical stability
    s0 = a*f(b)*f(c) / (epsilon + (f(a)-f(b))*(f(a)-f(c)))
    s1 = b*f(a)*f(c) / (epsilon + (f(b)-f(a))*(f(b)-f(c)))
    s2 = c*f(a)*f(b) / (epsilon + (f(c)-f(a))*(f(c)-f(b)))
    return s0+s1+s2

def s_secant(a, b, c):
    return b - (f(b)*(b-a)/(f(b)-f(a)))

def optimize():
    #define interval
    a = -.5 #if I let a<-1, it gives the wrong answer...let's just say brent's method doesn't work for
    #an interval with left bound < -1.
    b = .5 #b can't be too large either?
    tol = 1e-7
    #usually the first stage is to bracket. But our function f(x) is always > 0
    #so we can't use the usual f(a)*f(b) <= 0 condition
    if abs(f(a)) < abs(f(b)):
        a, b = b, a #swap bounds
    c = a
    flag = True #will keep track of whether previous step was bisection or not
    err = abs(b-a)
    err_list, b_list = [err], [b]
    while err > tol:
        if (f(a) != f(c)) or (f(b) != f(c)): #a,b,c are not collinear (pretty sure this is where
            # wikipedia gets it wrong?)
            s = s_quad_interp(a,b,c)
        else:
            #in this case, use secant method
            s = s_secant(a,b,c)
        #check for conditions (step size and tolerance checks)
        #if current step size > half the previous step size or if it is > half
        #the step size before the previous step (not necessary?), then do
        #bisection
        if ((flag == True) and (abs(s-b) >= abs(b-c)/2))\
           or ((flag == False) and (abs(s-b) >= abs(c-d)/2))\
           or ((flag == True) and (abs(b-c) < tol))\
           or ((flag == False) and (abs(c-d) < tol)):
            s = (a+b)/2 #bisection (can do golden ratio search instead)
            flag = True
        else:
            flag = False
        c, d = b, c
        if f(a)*f(s) < 0:
            b = s
        else:
            a = s
        if abs(f(a)) < abs(f(b)):
            a, b = b, a #swap if needed
        err = abs(b-a) #update error to check for convergence
        err_list.append(err)
        b_list.append(b)
    print(f'root = {b}')
    return b_list, err_list

def plot(b_list, err_list):
    log_err = [np.log10(err) for err in err_list]
    fig, axs = plt.subplots(2,1, sharex=True)
    ax0, ax1 = axs[0], axs[1]
    #plot root
    ax0.scatter(range(len(b_list)), b_list, marker = 'o', facecolor = 'red', edgecolor = 'k')
    ax0.plot(range(len(b_list)), b_list, 'r-', alpha = .5)
    ax1.plot(range(len(err_list)), log_err,'.-')
    ax1.set_xlabel('number of iterations')
    ax0.set_ylabel(r'$x_{min}$')
    ax1.set_ylabel(r'$\log{\delta}$')
    plt.savefig('convergence.png')
    
if __name__ == "__main__":
    b_list, err_list = optimize()
    plot(b_list, err_list)
