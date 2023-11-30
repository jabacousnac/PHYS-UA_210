from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy import linalg
import timeit

def read_file(filename):
    hdu_list = fits.open('specgrid.fits')
    logwave = hdu_list['LOGWAVE'].data
    flux = hdu_list['FLUX'].data
    return logwave, flux

def plot_galaxies(n, logwave, flux):
    fig, ax = plt.subplots()
    ax.set_xlabel(r'log$_{10}(\lambda [\AA])$')
    ax.set_ylabel(r'flux [erg s$^{-1}$ cm$^-2$ $\AA^{-1}$]')
    for i in range(n):
        ax.plot(logwave, flux[i], alpha = .5)
    plt.savefig('wave_v_spectrum.png', dpi=96) #should observe peak roughly @660 nm (red emission from n=3 to n=2)
    
def prep_flux(flux):
    '''preprocessing'''
    row_sum = flux.sum(axis=1)
    normed_flux = flux/row_sum[:, np.newaxis]
    row_ave = np.mean(normed_flux, axis = 1)
    normed_centered_flux = normed_flux - row_ave[:, np.newaxis]
    return normed_centered_flux, row_sum, row_ave

def PCA(logwave, flux):
    X = logwave
    Y, row_sum, row_ave = prep_flux(flux) #Y is normed_centered_flux, X is logwave

    #now, we want to construct the covariance matrix
    start01 = timeit.default_timer()
    C = np.transpose(Y) @ Y #should have shape NwavexNwave
    eigval, eigvec = linalg.eig(C)
    dt01 = timeit.default_timer()-start01
    n = 5
    fig, ax = plt.subplots()
    for i in range(n):
        ax.plot(logwave, np.real(eigvec[:, i]), alpha = .5)
    ax.set_xlabel(r'log$_{10}(\lambda [\AA])$')
    ax.set_ylabel('eigenvectors [arbs]')
    plt.savefig('eigenvectors_method01.png', dpi=96)

    #now do the SVD implementation, and then compare answers
    start02 = timeit.default_timer()
    U, W, V_T = linalg.svd(Y, full_matrices=False)
    W = np.diag(W)
    C_svd = np.transpose(V_T) @ W @ np.transpose(U) @ U @ W @ V_T
    eigval, eigvec = linalg.eig(C_svd)
    dt02 = timeit.default_timer()-start02
    fig, ax = plt.subplots()
    for i in range(5):
        ax.plot(logwave, np.real(eigvec[:, i]), alpha = .5)
    ax.set_xlabel(r'log$_{10}(\lambda [\AA])$')
    ax.set_ylabel('eigenvectors [arbs]')        
    plt.savefig('eigenvectors_method02.png', dpi=96)
    print('method01: {}s'.format(dt01) + '\n' + 'method02: {}s'.format(dt02))

    #find the condition number (use L-2 norm)
    norm_C = linalg.norm(C, 2)
    norm_Cinv = linalg.norm(linalg.inv(C))
    print(f'kappa_C = {norm_C*norm_Cinv}') #see how kappa_C is extremely high, which means it is
    #almost singular. Thus, SVD would be a better way to construct the covariance matrix because
    #of numerical stability!

    #find the condition number of R matrix
    norm_R = linalg.norm(np.diag(Y), 2)
    norm_Rinv = linalg.norm(1./np.diag(Y),2)
    print(f'kappa_R = {norm_R*norm_Rinv}') #still high, but less
    #I get kappa_C = 16912535552, kappa_R = 420990

    #now we want to reconstruct the spectrum
    _, eigvec = linalg.eig(C_svd)
    #take the first l eigenvectors (l=5)
    v = np.transpose(eigvec[:5]) #4001x5
    c_i = Y @ v #9713x5
    c_0 = c_i[:,0]
    c_1 = c_i[:,1]
    c_2 = c_i[:,2]
    #plot the c_i's
    fig, axs = plt.subplots(2,1,sharex=True)
    axs[0].scatter(c_0, c_1, color = 'red', alpha = .5)
    axs[1].scatter(c_0, c_2, color = 'red', alpha = .5)
    axs[1].set_xlabel(r'$c_{0}$')
    axs[0].set_ylabel(r'$c_{1}$')
    axs[1].set_ylabel(r'$c_{2}$')
    plt.tight_layout()
    plt.savefig('principal_components.png', dpi = 96)

    #print some shapes so we know how to reconstruct the spectrum
    #print('shapes:', c_i.shape, v.shape, row_sum.shape, row_ave.shape, 'shapes_end')
    dim = v.shape[0]
    sums = np.array([dim*[q] for q in row_sum])
    aves = np.array([dim*[q] for q in row_ave])
    reconstructed = (c_i @ v.T) + aves
    print(reconstructed.shape)

    #for the last part, I got a little lazy. Loop over a range of 20, and each time
    #reconstruct the OG dataset by dotting the coeffs and the eigenvectors. Compute
    #error by comparing with the OG dataset. 

if __name__ == "__main__":
    logwave, flux = read_file('specgrid.fits')
    #plot_galaxies(3, logwave, flux)
    PCA(logwave, flux)
