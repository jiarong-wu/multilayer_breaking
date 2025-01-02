""" Utility functions for spectrum manipulation. """
import numpy as np
from scipy.interpolate import griddata

""" General conversion between cartesian and polar coordinate. """
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)


def spectrum_integration_log(eta, L, N, CHECK=False, method='nearest'):   
    """ This function performs azimuthal integration of the 2D spectrum.
        When in doubt, enable CHECK so that the integration is printed out at each step to make sure that 
        units are consistent and we always recover the variance of the eta matrix (wave height). 
        The output k is log-spaced instead of linearly-spaced as its predecessor. """  
    
    if CHECK: print (np.var(eta))
    spectrum = np.fft.fft2(eta) / (N*N)**0.5 # FFT normalization 
    F = np.absolute(spectrum)**2 / N**2 # Per area normalization
    if CHECK: print (np.sum(F))
    wavenumber = 2*np.pi*np.fft.fftfreq(n=N,d=L/N)
    kx = np.fft.fftshift(wavenumber); ky = kx
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    
    ''' Setting the wavenumber interpolated to '''
    theta = np.linspace(-N/2, N/2, 361)/(N)*2*np.pi
    k = np.logspace(-2, np.log10(1.414*kx.max()), 200) # get resolution higher, and in logspace
    dk = np.roll(k, -1) - k # now dk is varying
    dk = np.array(dk); dk[-1]=0    
    dtheta = theta[1]-theta[0]   
    dkx = kx[1] - kx[0]; dky = ky[1] - ky[0]

    F_center = np.fft.fftshift(F)/dkx/dky # Further normalization by independent variables
    k_tile, theta_tile = np.meshgrid(k,theta)
    kxp_tile, kyp_tile = pol2cart(k_tile, theta_tile)
    
    if method == 'linear':
        F_center_polar = griddata((kx_tile.ravel(),ky_tile.ravel()), F_center.ravel(),
                                  (kxp_tile, kyp_tile), method='linear', fill_value=0)
    elif method == 'nearest':
        F_center_polar = griddata((kx_tile.ravel(),ky_tile.ravel()), F_center.ravel(), 
                              (kxp_tile, kyp_tile), method='nearest') # It's faster if using nearest, and seems as accurate
        
    F_center_polar_integrated = np.sum(F_center_polar*k_tile, axis=0)*dtheta # Azimuthal integration
    if CHECK: print (np.sum(F_center_polar_integrated*dk))
        
    return k, theta, F_center_polar_integrated


def spectrum_integration_linear(eta, CHECK=False, N=512, L=50):
    """ This function performs azimuthal integration of the 2D spectrum.
        When in doubt, enable CHECK so that the integration is printed out at each step to make sure that 
        units are consistent and we always recover the variance of the data. 
        This is the old linearly spaced version. """  
    
    if CHECK: print (np.var(eta))
    spectrum = np.fft.fft2(eta) / (N*N)**0.5 # FFT normalization 
    F = np.absolute(spectrum)**2 / N**2 # Per area normalization
    if CHECK: print (np.sum(F))
    wavenumber = 2*np.pi*np.fft.fftfreq(n=N,d=L/N)
    kx = np.fft.fftshift(wavenumber); ky = kx
    kx_tile, ky_tile = np.meshgrid(kx,ky)
    
    theta = np.arange(-N/2,N/2)/(N)*2*np.pi   
    k = wavenumber[0:int(N/2)]
    dkx = kx[1] - kx[0]; dky = ky[1] - ky[0]
    dk = k[1]-k[0]; dtheta = theta[1]-theta[0]
    
    F_center = np.fft.fftshift(F)/dkx/dky # Further normalization by independent variables
    k_tile, theta_tile = np.meshgrid(k,theta)
    kxp_tile, kyp_tile = pol2cart(k_tile, theta_tile)
    F_center_polar = griddata((kx_tile.ravel(),ky_tile.ravel()), F_center.ravel(), (kxp_tile, kyp_tile), method='linear', fill_value=0)
    F_center_polar_integrated = np.sum(F_center_polar*k_tile, axis=0)*dtheta # Azimuthal integration
    if CHECK: print (np.sum(F_center_polar_integrated)*dk)
    return k, theta, F_center_polar_integrated

def steepness_trunc_non_uniform (F,k):
    """ This function computes the cumulative steepness mu and cumulative wave height Hs from
        the input uni-directional F and k. k can be non-linearly spaced. Output mu and Hs are of the same
        dimension as input k. """
    mu = np.zeros(len(k))
    Hs = np.zeros(len(k))
    dk = np.roll(k, -1) - k
    dk = np.array(dk); dk[-1]=0
    for i,k_ in enumerate(k):
        mu[i] = np.sum(k[:i]**2*F[:i]*dk[:i])**0.5
        Hs[i] = 4*np.sum(F[:i]*dk[:i])**0.5
    return mu, Hs

def eta_random(kx_tile, ky_tile, F_kxky_tile, dkx, dky, omega_tile, x_tile, y_tile, t=0):
    """ This function generates a random-phased eta filed from a given kx-ky spectrum.
        Arguments:
            kx_tile, ky_tile, F_kxky_tile, dkx, dky, omega_tile: same dimensions 
            x_tile, y_tile: discretized location x and y, can be generated using this example code
                            x = np.linspace(-L/2,L/2,N_grid,endpoint=False)+L/N_grid/2 
                            y = np.linspace(-L/2,L/2,N_grid,endpoint=False)+L/N_grid/2
                            x_tile, y_tile = np.meshgrid(x, y)          
            t: the time that the eta field is computed, by default zero       
    """
    np.random.seed(0) 
    phase_tile = np.random.random_sample(kx_tile.shape)*2*np.pi # Add a random phase field
    eta_tile = np.zeros(x_tile.shape)
    for i1 in range(0,kx_tile.shape[0]):
        for i2 in range(0,kx_tile.shape[1]):
            ampl = (2*F_kxky_tile*dkx*dky)**0.5
            # How to exactly represent integrate over dk_x*dk_y*eta_hat?
            # mode = (F_kdirectional**0.5)*np.cos((kx_tile*x_tile[i1,i2]+ky_tile*y_tile[i1,i2])+
            #                                    phase_tile)*(kmod_tile*(kmod[1]-kmod[0])*(theta[1]-theta[0])) # uniform space in k and theta
            a = (kx_tile*x_tile[i1,i2]+ky_tile*y_tile[i1,i2])-omega_tile*t+phase_tile
            mode = ampl*(np.cos(a)) # uniform spacing in kx and ky
            eta_tile[i1,i2] = np.sum(mode)    
    return eta_tile, phase_tile

def eta_focusing(kx_tile, ky_tile, F_kxky_tile, dkx, dky, omega_title, x_tile, y_tile, t=0, tb = 40, xb = 0, yb = 0):
    """ This function generates a eta filed from a given kx-ky spectrum, and the phase is initialized so that all the modes
        linearly focus at location (xb, yb) at time tb.
        Arguments are the same with eta_random. Extra arguments are tb, xb, yb, which control the focusing.    
    """
    phase_tile = -kx_tile*xb-ky_tile*yb+omega_tile*tb
    eta_tile = np.zeros(x_tile.shape)
    for i1 in range(0,kx_tile.shape[0]):
        for i2 in range(0,kx_tile.shape[1]):
            ampl = (2*F_kxky_tile*dkx*dky)**0.5
            a = (kx_tile*x_tile[i1,i2]+ky_tile*y_tile[i1,i2])-omega_tile*t+phase_tile
            mode = ampl*(np.cos(a)) # uniform space in kx and ky
            eta_tile[i1,i2] = np.sum(mode)  
    return eta_tile, phase_tile