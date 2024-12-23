from scipy.interpolate import interp1d
import xarray as xr
import xgcm
import gc
import numpy as np

class interp1d_class():
    def __init__(self, znew, fill_value):
        self.znew = znew
        self.fill_value = fill_value
    def __call__(self, x, y):
        f = interp1d(x, y, kind='linear', fill_value=self.fill_value, bounds_error=False)
        return f(self.znew)
    
''' fill_value: array-like or (array-like, array_like) or 'extrapolate' '''
interpz = lambda znew, z, ds, fill_value: xr.apply_ufunc (
    interp1d_class(znew, fill_value),
    #lambda x, y: interpolate_1d(x, y, znew),  # The function to apply if interpolate_1d is not defined as a class
    z,  # Original vertical coordinates
    ds,  # Data (Dataarray or Dataset) to interpolate
    input_core_dims=[['zl'], ['zl']],  # Core dimensions for each input
    output_core_dims=[['zl']],  # Core dimensions for the output
    exclude_dims=set(('zl',)),
    output_sizes={'zl':len(znew)},
    output_dtypes=['float32'],
    vectorize=True,  # Enable vectorization
    dask="parallelized",  # Parallelize using Dask if the data is large
)

def dissipation_layer(ds, grid):
    
    delta = ds.x[1]-ds.x[0]
    
    dudx = grid.interp(grid.diff(ds.ux, 'X'), 'X')/delta
    dudy = grid.interp(grid.diff(ds.ux, 'Y'), 'Y')/delta
    dudzl = grid.interp(grid.diff(ds.ux, 'Z'), 'Z')
    dvdx = grid.interp(grid.diff(ds.uy, 'X'), 'X')/delta
    dvdy = grid.interp(grid.diff(ds.uy, 'Y'), 'Y')/delta
    dvdzl = grid.interp(grid.diff(ds.uy, 'Z'), 'Z')
    dwdx = grid.interp(grid.diff(ds.uz, 'X'), 'X')/delta
    dwdy = grid.interp(grid.diff(ds.uz, 'Y'), 'Y')/delta
    dwdzl = grid.interp(grid.diff(ds.uz, 'Z'), 'Z')
    
    dzdx = grid.interp(grid.diff(ds.z, 'X'), 'X')/delta
    dzdy = grid.interp(grid.diff(ds.z, 'Y'), 'Y')/delta
    dzdzl = grid.interp(grid.diff(ds.z, 'Z'), 'Z')
    
    ds['dudz'] = dudzl/dzdzl
    ds['dudy'] = dudy - ds['dudz']*dzdy
    ds['dudx'] = dudx - ds['dudz']*dzdx
    ds['dvdz'] = dvdzl/dzdzl
    ds['dvdy'] = dvdy - ds['dvdz']*dzdy
    ds['dvdx'] = dvdx - ds['dvdz']*dzdx
    ds['dwdz'] = dwdzl/dzdzl
    ds['dwdy'] = dwdy - ds['dwdz']*dzdy
    ds['dwdx'] = dwdx - ds['dwdz']*dzdx
    
    del(dudx, dudy, dudzl, dvdx, dvdy, dvdzl, dwdx, dwdy, dwdzl, dzdx, dzdy, dzdzl)
    gc.collect()
    
    ds['epsilon'] = 2*(ds.dudx**2 + 2*((ds.dudy + ds.dvdx)/2.)**2 + 2*((ds.dudz+ds.dwdx)/2.)**2 + ds.dvdy**2 + 2*((ds.dvdz+ds.dwdy)/2.)**2 + ds.dwdz**2)
    
    ds['omegaxp'] = ds.dwdy - ds.dvdz
    ds['omegayp'] = ds.dudz - ds.dwdx
    ds['omegazp'] = ds.dvdx - ds.dudy
    ds['vort2'] = ds.omegaxp**2 + ds.omegayp**2 + ds.omegazp**2
    
    return ds

def compute_diss (ds):
    grid = xgcm.Grid(ds, 
                 coords={
                     'X':{'center':'x', 'left':'x_g'},
                     'Y':{'center':'y', 'left':'y_g'},
                     'Z':{'center':'zl', 'left':'zl_g'},
                 },
                 periodic={'X':'True','Y':'True','Z':'False'},
                 boundary={'Z':'fill'},
                 fill_value={'Z':0})
    ds = dissipation_layer(ds, grid)
    return ds

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