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

##### Converting different mesh format #####

""" Assemble the 3D array of vertice for pcolormesh. We construct the x_mesh, y_mesh, and z_mesh 
    (of dimension Nl+1, Nh+1, Nh+1) based on h field (and the horizontal dimensions L0), 
    while the other field remains unchanged. """
def array_to_mesh (h_ensem, L0=200, H=40, Nh=512, Nl=15):
 
    x_mesh = np.zeros([Nl+1,Nh+1,Nh+1]) # Different from the vtk format
    y_mesh = np.zeros([Nl+1,Nh+1,Nh+1])
    z_mesh = np.zeros([Nl+1,Nh+1,Nh+1])

    h_ensem_expand = np.zeros([Nl,Nh+1,Nh+1]) # Need to go from centered to grid, pad the array
    h_ensem_expand[:,:Nh,:Nh] = np.copy(h_ensem) # Need to go from centered to grid
    h_ensem_expand[:,Nh,:Nh] = np.copy(h_ensem[:,Nh-1,:Nh])
    h_ensem_expand[:,:Nh,Nh] = np.copy(h_ensem[:,:Nh,Nh-1])
    h_ensem_expand[:,Nh,Nh] = np.copy(h_ensem[:,Nh-1,Nh-1])
    h_ensem_expand = np.array(h_ensem_expand)

    xarray = np.linspace(-L0/2, L0/2, Nh+1, endpoint=True)
    yarray = np.linspace(-L0/2, L0/2, Nh+1, endpoint=True)
    
    for k in range(Nl+1):
        for i in range(Nh+1):
            for j in range(Nh+1):
                z_mesh[k,i,j] = np.sum(h_ensem_expand[:k,i,j]) - H
                x_mesh[k,i,j] = xarray[i]
                y_mesh[k,i,j] = yarray[j]

    return x_mesh, y_mesh, z_mesh


""" convert the 3D array to vtk file for paraview. Need to specify L0 and H """
from pyevtk.hl import gridToVTK
def array_to_vtk (h_ensem, ux_ensem, uy_ensem, uz_ensem, ichoice, 
                  filepath='/projects/DEIKE/jiarongw/multilayer/paraview/vtk/', 
                  L0=200, H=40, Nh=512, Nl=15):
    fieldnames = ['x', 'y', 'z', 'ux', 'uy', 'uz', 'f']

    x_vtk = np.zeros([Nh+1,Nh+1,Nl+1])
    y_vtk = np.zeros([Nh+1,Nh+1,Nl+1])
    z_vtk = np.zeros([Nh+1,Nh+1,Nl+1])
    f_vtk = np.zeros([Nh,Nh,Nl])
    ux_vtk = np.zeros([Nh,Nh,Nl])
    uy_vtk = np.zeros([Nh,Nh,Nl])
    uz_vtk = np.zeros([Nh,Nh,Nl])

    h_ensem_expand = np.zeros([Nl,Nh+1,Nh+1]) # Need to go from centered to grid, pad the array
    h_ensem_expand[:,:Nh,:Nh] = np.copy(h_ensem) # Need to go from centered to grid
    h_ensem_expand[:,Nh,:Nh] = np.copy(h_ensem[:,Nh-1,:Nh])
    h_ensem_expand[:,:Nh,Nh] = np.copy(h_ensem[:,:Nh,Nh-1])
    h_ensem_expand[:,Nh,Nh] = np.copy(h_ensem[:,Nh-1,Nh-1])
    h_ensem_expand = np.array(h_ensem_expand)

    xarray = np.linspace(-L0/2, L0/2, Nh+1, endpoint=True)
    yarray = np.linspace(-L0/2, L0/2, Nh+1, endpoint=True)

    for k in range(Nl+1):
        for i in range(Nh+1):
            for j in range(Nh+1):
                z_vtk[i,j,k] = np.sum(h_ensem_expand[:k,i,j]) - H
                x_vtk[i,j,k] = xarray[i]
                y_vtk[i,j,k] = yarray[j]

    for k in range(Nl):
        for i in range(Nh):
            for j in range(Nh):
                ux_vtk[i,j,k] = ux_ensem[k,i,j]
                uy_vtk[i,j,k] = uy_ensem[k,i,j]
                uz_vtk[i,j,k] = uz_ensem[k,i,j]
                if k == Nl-1: # surface layer
                    f_vtk[i,j,k] = 0
                else:
                    f_vtk[i,j,k] = 1
                    
    gridToVTK(filepath + "structured_%g" %ichoice, x_vtk, y_vtk, z_vtk, cellData = {"f": f_vtk, "ux": ux_vtk, "uy": uy_vtk, "uz": uz_vtk})
    return x_vtk, y_vtk, z_vtk