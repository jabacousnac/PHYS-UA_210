import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as lg
import pandas as pd
import scipy.signal as sg
import scipy.optimize as opt

def plotter(X, Y):
    fig, ax = plt.subplots()
    ax.grid(alpha=.5)
    ax.scatter(X, Y, color = 'r', alpha = .5)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('signal [arbs]')
    plt.savefig('signal_v_time.png', dpi=96)
    
def get_data(filename):
    df = pd.read_csv(filename, delimiter='|', header=0, names=['a', 'time', 'signal', 'b']) #the file is formatted such that the two delimeters on either side mean that there will be four columns
    X = df['time']
    Y = df['signal']
    return X,Y

def poly_fitter(X, Y, n):
    X = np.array(X).flatten()
    Y = np.array(Y).flatten()
    A = np.vander(X, n+1) #create the design matrix
    U, W, V_T = lg.svd(A, full_matrices=False)
    #print(U, W, V_T)
    W_inv = np.diag(1./W) #W is diagonal -> it's inverse should be 1/W
    coeffs = np.transpose(V_T) @ W_inv @ np.transpose(U) @ Y
    #since V_T@V = 1 = U_T@U, we can multiply both sides, to get out the coeffs as above^
    c3, c2, c1, c0 = coeffs
    X_fit = np.linspace(min(X), max(X), 100)
    Y_fit = [c0 + c1*x + c2*x**2 + c3*x**3 for x in X_fit]
    Y_hat = [c0 + c1*x + c2*x**2 + c3*x**3 for x in X.flatten()]
    #calculate residuals:
    res = np.array(Y)-np.array(Y_hat)
    sres = np.sum(res**2)
    stot = np.sum((np.array(Y)-np.mean(Y))**2)
    print('r_sq = {}'.format(1-sres/stot)) #rsq
    #plot:
    fig, ax = plt.subplots()
    ax.grid(alpha = .5)
    ax.scatter(X, Y, c='blue', alpha = .5, label='data')
    ax.plot(X_fit, Y_fit, 'k--', alpha = 1., label='poly (n={}) fit\n'.format(n) + r'$r^2={}$'.format(np.round(1-sres/stot,2)))
    ax.set_xlabel('time [s]')
    ax.set_ylabel('signal [arbs.]')
    ax.legend()
    plt.tight_layout()
    plt.savefig('signal_v_time_poly_fit_n={}.png'.format(n), dpi=96)

def lomb_scargle(X, Y):
    #do lomb-scargle periodgram
    X = np.array(X)
    Y = np.array(Y)
    w = np.linspace(min(X), max(X), 50)
    pgram = sg.lombscargle(X.flatten(), Y.flatten(), w, normalize=False)
    fig, ax = plt.subplots()
    ax.plot(w, pgram, c='orange')
    ax.set_xlabel(r'$\omega$ [rad/s]')
    ax.set_ylabel(r'Amplitude [arbs.]')
    plt.savefig('periodgram.png', dpi=96)

def get_trig_coeffs(n, X, Y, P): #P is the period
    X = np.array(X).flatten()
    Y = np.array(Y).flatten()
    A = np.zeros((len(X), n+1))
    A[:, -1] = 1
    for q in range(0, n):
        A[:, q] = np.sin(2*np.pi*q*X/P) #only doing sine, but you can likewise
        #expand decision matrix to include cosine
    U, W, V_T = lg.svd(A, full_matrices=False)
    idx = np.where(W<1e-16)
    W[idx] = 1. #ehh
    W_inv = np.diag(1./W)
    coeffs = np.transpose(V_T) @ W_inv @ np.transpose(U) @ Y
    return coeffs

def trig_fitter(X, Y):
    X = np.array(X).flatten()
    Y = np.array(Y).flatten()
    X_fit = np.linspace(min(X), max(X), 100)
    P = max(X)/2
    pmin = max(X)/10
    #start plotting:
    fig, ax = plt.subplots()
    ax.scatter(X, Y, color='red', alpha = .2)
    ax.set_xlabel('time [s]')
    ax.set_ylabel('signal [arbs.]')
    n=5
    while P > pmin:
        a, b, c, d, e, f = get_trig_coeffs(n, X, Y, P)
        print(a,b,c,d,e,f)
        Y_fit = [a*np.sin(n*np.pi*x/P) + b*np.sin((n-1)* np.pi*x/P) + c*np.sin((n-2)*np.pi*x/P) + d*np.sin((n-3)*np.pi*x/P) + e*np.sin(np.pi*x/P) + f for x in X_fit]
        Y_hat = [a*np.sin(n*np.pi*x/P) + b*np.sin((n-1)* np.pi*x/P) + c*np.sin((n-2)*np.pi*x/P) + d*np.sin((n-3)*np.pi*x/P) + e*np.sin(np.pi*x/P) + f for x in X]        
        ax.plot(X_fit, Y_fit, label=r'$P/P_0 = {}$'.format(np.around(P/max(X),4)))
        #calculate residuals:
        res = Y-np.array(Y_hat)
        sres = np.sum(res**2)
        stot = np.sum((Y-np.mean(Y))**2)
        print('r_sq = {}'.format(1-sres/stot)) #rsq
        P-=pmin
    ax.legend()
    plt.tight_layout()
    plt.savefig('trig_fitting.png', dpi=96)
    #disclaimer: the trig portion of my solution isnt quite right. You'd need a cosine fourier series should the function be even. Should also go to higher harmonic. Sorry, it's been a busy week :(
    
if __name__ == '__main__':
    out = get_data('signal.dat')
    #plotter(out[0], out[1])
    #poly_fitter(out[0], out[1], 10)
    trig_fitter(out[0], out[1])
