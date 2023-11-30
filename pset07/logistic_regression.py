import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize

def nLL(beta, X_data, Y_data):
    #return the negative log likelihood, which we can then minimize
    '''
    log-likelihood becomes the sum:
    i = 0 to N (where N is total amount of datapoints)
    Sum (y_i lnP + (1-y_i)ln(1-P))
    ^https://arunaddagatla.medium.com/maximum-likelihood-estimation-in-logistic-regression-f86ff1627b67
    '''
    # probability distribution is logistic function:
    p = lambda x, beta: 1 / ((1 + np.exp(-beta[0] + x * beta[1])))
    p_i = p(X_data, beta)
    epsilon = 1e-8 #for numerical stability
    ll_sum = np.sum(Y_data@np.log(p_i+epsilon) + (1-Y_data)@np.log(1.-p_i+epsilon))
    return -ll_sum

def minimize_nLL(X, Y):
    x0 = np.array([-5, 0])
    res = minimize(nLL, x0=x0, args = (X, Y), method='BFGS')
    return res

def plot_res(X, Y):
    #minimize neg log likelihood
    res = minimize_nLL(X, Y)
    beta0_hat, beta0_hat_err = res.x[0], np.sqrt(res.hess_inv[0,0])
    beta1_hat, beta1_hat_err = res.x[1], np.sqrt(res.hess_inv[1,1])
    #plot data
    fig, ax = plt.subplots()
    ax.scatter(X, Y, marker='o', facecolor='r', edgecolor='k')
    ax.set_xlabel('age')
    ax.set_ylabel('ans')
    #plot fit
    p = lambda x, beta: 1 / ((1 + np.exp(-beta[0] + x * beta[1])))
    X_fit = np.linspace(min(X), max(X), 100)
    Y_fit = p(X_fit, [beta0_hat, beta1_hat])
    params = r'$\beta_0 $ = {} +/-'.format(np.around(beta0_hat,3)) + f'{beta0_hat_err:,.2f}' + '\n'+\
             r'$\beta_1 $ = {} +/-'.format(np.around(beta1_hat,3)) + f'{beta1_hat_err:,.2f}'
    ax.plot(X_fit, Y_fit, 'b--', label=params)
    ax.legend()
    plt.savefig('logistic_curve_fit.png', dpi=96)

if __name__ == "__main__":
    #load data
    df = pd.read_csv('survey.csv')
    X, Y = np.array(np.array(df['age'])), np.array(df['recognized_it'])
    plot_res(X, Y)

    
