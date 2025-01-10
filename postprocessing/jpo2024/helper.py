import numpy as np
import pandas as pd
import os
from mlpython.specfunc import get_mss_Hs_spectrum

''' Figure 1: compute time binned dEdt and wave characteristics '''
def compute_dEdt(ds, tbins):
    dEkdt = []; dEpdt = []            
    for i in range(len(tbins[:-1])):
        dEk = ds.ke.sel(t=tbins[i+1], method='nearest') - ds.ke.sel(t=tbins[i], method='nearest')    
        dEp = ds.gpe.sel(t=tbins[i+1], method='nearest') - ds.gpe.sel(t=tbins[i], method='nearest')
        dEkdt.append(dEk/(tbins[i+1]-tbins[i])); dEpdt.append(dEp/(tbins[i+1]-tbins[i]))
    dEkdt = np.array(dEkdt); dEpdt = np.array(dEpdt)
    dEdt = (dEkdt + dEpdt)/ds.attrs['L']**2
    return dEdt

def compute_waves(ds, tbins):
    mus = []; Hss = []
    #### Use the middle of the bin ####
    for i in range(len(tbins[:-1])):
        tc = (tbins[i] + tbins[i+1])/2
        mu, Hs, kmod, Fkmod = get_mss_Hs_spectrum (ds.eta.sel(t=tc, method='nearest').values, 
                                                   ds.attrs['L'], 2**ds.attrs['LEVEL'])   
        mus.append(mu); Hss.append(Hs)  
    return mus, Hss


''' For figure 4 plot literature data '''
def plot_literature(ax, path):
    # Sutherland 2015
    for i in range(1,40):
        file = path + 'literature_data/SM15/SoCal2010_C1_B2_%d.txt' %i
        if os.path.isfile(file):
            S15 = pd.read_csv(file, delimiter=' ')
            ax.plot(S15.epsilon, S15.z, '*', c='gray')
        file = path + 'literature_data/SM15/SoCal2010_C1_B3_%d.txt' %i
        if os.path.isfile(file):
            S15 = pd.read_csv(file, delimiter=' ')
            ax.plot(S15.epsilon, S15.z, '*', c='gray')    
        file = path + 'literature_data/SM15/SoCal2010_C1_VERT_%d.txt' %i
        if os.path.isfile(file):
            S15 = pd.read_csv(file, delimiter=' ')
            ax.plot(S15.epsilon, S15.z, '-*', c='gray') 
            
    for i in range(4,10):
        S15 = pd.read_csv(path + 'literature_data/SM15/SoCal2010_C2_B2_%d.txt' %i, delimiter=' ')
        ax.plot(S15.epsilon, S15.z, '*', c='gray')   
    for i in (2,4,5,6,7):    
        S15 = pd.read_csv(path + 'literature_data/SM15/SoCal2010_C3_B2_%d.txt' %i, delimiter=' ')
        ax.plot(S15.epsilon, S15.z, '*', c='gray')   

    S15 = pd.read_csv(path + 'literature_data/SM15/SoCal2010_C3_B2_%d.txt' %7, delimiter=' ')
    ax.plot(S15.epsilon, S15.z, '*', c='gray', label='SM15')   

    # Other data
    SL03 = pd.read_csv(path + 'literature_data/SL03.csv', names=['x', 'y'])
    ax.plot(SL03.x, SL03.y, '-', c='k', label='SL03')

    T96fit = pd.read_csv(path + 'literature_data/T96fit.csv', names=['x', 'y'])
    # axes[1].plot(T96fit.x, T96fit.y, '--', c='k', label='T96 fit')
    # Found the expression for the fit
    ax.plot(T96fit.y**(-2)*0.3, T96fit.y, '--', c='k', label='T96 fit')
    ax.vlines(x=0.83, ymin=-0.6, ymax=-0.1, linestyle='--', color='k')

    ax.plot(12.0484, -0.2830, 'o', c='k')
    T96 = pd.read_csv(path + 'literature_data/T96.csv', names=['x', 'y'])
    ax.plot(T96.x, T96.y, 'o', c='k', label='T96')

    xstart = 0.01; ystart = -2.
    ax.plot(T96fit.x[0:8], -T96fit.x[0:8]**(-0.5)/(-T96fit.x[0]**(-0.5))*(ystart), c='gray', alpha=0.8)
    ax.annotate('$z^{-2}$',(T96fit.x[8],-T96fit.x[8]**(-0.5)/(-T96fit.x[0]**(-0.5))*(ystart)), fontsize=6)
    ax.plot(T96fit.x[0:8], -T96fit.x[0:8]**(-1)/(-T96fit.x[0]**(-1))*(ystart), c='gray', alpha=0.8)
    ax.annotate('$z^{-1}$',(T96fit.x[8],-T96fit.x[8]**(-1)/(-T96fit.x[0]**(-1))*(ystart)), fontsize=6)

    D96young = pd.read_csv(path + 'literature_data/D96young.csv', names=['x', 'y'])
    ax.plot(D96young.x, D96young.y, '^', c='k', label='D96')
    D96old = pd.read_csv(path + 'literature_data/D96old.csv', names=['x', 'y'])
    ax.plot(D96old.x, D96old.y, '^', c='k')

    AM95 = pd.read_csv(path + 'literature_data/AM95.csv', names=['x', 'y'])
    ax.plot(AM95.x, AM95.y, 's', c='k', label='AM95')

    # S15 = pd.read_csv(path + 'literature_data/Sutherland_upper.txt', names=['x', 'y'])
    # axes[1].plot(S15.x, S15.y, '*', c='k', label='S15')
    # S15 = pd.read_csv(path + 'literature_data/Sutherland_lower.txt', names=['x', 'y'])
    # axes[1].plot(S15.x, S15.y, '*', c='k')

    # After getting data from Peter Sutherland
    # for i in range(1,5):
    #     S15 = pd.read_csv(path + 'literature_data/HIRES2010_%d.txt' %i, delimiter=' ')
    #     plt.plot(S15.epsilon, S15.z, '*', c='k', label='S15')
    # for i in (10,11,13,14):
    #     S15 = pd.read_csv(path + 'literature_data/RaDyO2019_%d.txt' %i, delimiter=' ')
    #     plt.plot(S15.epsilon, S15.z, '*', c='k', label='S15')
    
''' Figure 4 axis locator? '''
from matplotlib.ticker import Locator
class MinorSymLogLocator(Locator):
    """
    Dynamically find minor tick positions based on the positions of
    major ticks for a symlog scaling.
    """
    def __init__(self, linthresh):
        """
        Ticks will be placed between the major ticks.
        The placement is linear for x between -linthresh and linthresh,
        otherwise its logarithmically
        """
        self.linthresh = linthresh

    def __call__(self):
        'Return the locations of the ticks'
        majorlocs = self.axis.get_majorticklocs()

        # iterate through minor locs
        minorlocs = []

        # handle the lowest part
        for i in range(1, len(majorlocs)):
            majorstep = majorlocs[i] - majorlocs[i-1]
            if abs(majorlocs[i-1] + majorstep/2) < self.linthresh:
                ndivs = 10
            else:
                ndivs = 9
            minorstep = majorstep / ndivs
            locs = np.arange(majorlocs[i-1], majorlocs[i], minorstep)[1:]
            minorlocs.extend(locs)

        return self.raise_if_exceeds(np.array(minorlocs))

    def tick_values(self, vmin, vmax):
        raise NotImplementedError('Cannot get tick locations for a '
                                  '%s type.' % type(self))

''' For figure 4 compute color of lines based on sigma, given the time series of eta '''  
import xarray as xr
import matplotlib as mpl  
def compute_sigma_color (labels, tbin=[160,180], path='/Users/jiarongw/Data/multilayer_data/JPO2024/processed/'):
    sigmas = []
    ccs = []
    Hss = []
    for label in labels:
        filename = path + label + '/series.nc'
        ds = xr.open_dataset(filename, engine='h5netcdf')
        sigma, Hs = compute_waves(ds, tbin) # Compute sigma and Hs, with only one bin
        cc = mpl.colormaps['Oranges']((sigma[0]-0.05)/(0.2-0.05)) 
        sigmas.append(sigma[0]); Hss.append(Hs[0]); ccs.append(cc)        
    return sigmas, Hss, ccs

''' Figure 5 compute Phi '''
def compute_Phi (labels=['C1','C2','C3','C4','C5'], tcs=[110,130,150,170], 
                 basepath='/Users/jiarongw/Data/multilayer_data/JPO2024/processed/'):
    dicts = []
    for label in labels:
        phis = []
        filename = basepath + label + '/epsilon1d.nc'
        ds = xr.open_dataset(filename, engine='h5netcdf')
        for tc in tcs:
            phi = ds.sel(t=tc, method='nearest').epsilon.integrate('z')
            phis.append(phi.values)
        ds.close()
        phi_dict = {'label':label, 'tcs':tcs, 'phis':np.array(phis).flatten()}
        dicts.append(phi_dict)
    return dicts