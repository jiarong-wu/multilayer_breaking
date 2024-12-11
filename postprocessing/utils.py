import numpy as np

###### Compute Stokes drift given a spectrum #######

def spectrum_us(k, F, zarray):   
    """ Given the spectrum and the z coordinate, compute the stokes drift. """  
    
    us = np.zeros(len(zarray))
    dk = np.roll(k, -1) - k # now dk is varying
    dk = np.array(dk); dk[-1]=0    
    g = 9.8
    for i,z in enumerate(zarray):
        us[i] = 2*g**0.5*np.sum(k**1.5*F*dk*np.exp(2*k*z))
        
    return us