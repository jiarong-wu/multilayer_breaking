import xarray as xr
import numpy as np
import pickle
from mlpython.utils import array_to_mesh

''' Basic reading files '''
def read (filepath='/projects/DEIKE/jiarongw/multilayer/stokes/stokes_8_20_Htheta0.51/field_t5/',
          pre='eta_matrix_', index=0, N=512):
    filename = filepath + pre + '%g' %index
    f = np.fromfile(filename, dtype=np.float32)
    f = f.reshape(N+1,N+1); f = f[1:,1:]
    return f

def read_t(fieldnames=['h','ux','uy','uz'], t=5, Nh=30, Nl=512, path='/projects/DEIKE/jiarongw/multilayer/revision/stokes_ml_1/'):
    ''' Arguments:
            fieldnames: field names list
            Nh: horizontal resolution
            NL: layer numbers 
            path: path to the main folder        
        Returns:
            fields corresponding to fieldnames (Nl*Nx*Ny)
            Notice that in some cases the first layer was not written correctly and requires further filtering.
        '''
    folder = path + 'field/'
    fields = []
    for fieldname in fieldnames:
        field = []
        """ axis0-z; axis1-x; aixs2-y"""
        for l in range (0, Nl):
            fieldl = read(filepath = folder, pre = fieldname + '_matrix_t%g_l' %t, index = l, N = Nh)
            field.append(fieldl)
        fields.append(np.array(field))
    return fields


''' Function to read files and assemble them into an xarray dataset '''
def read_netcdf (field_path, config, t, fieldnames=['h','ux','uy','uz','omegax','omegay','omegaz','dzdx','dzdy','dzdxc','dzdyc']):
    # 3D field path and 2D eta path relative to the main path
#     field_path = path + 'vort/'
#     eta_path = path + 'surface/'
    # Reading fields (centered) according to the name list
    fields = read_t(fieldnames=fieldnames, t=t, Nh=2**config['LEVEL'], Nl=config['NL'], path=field_path)    
    # 3D meshes at vertices
    x_mesh, y_mesh, z_mesh = array_to_mesh (fields[fieldnames.index('h')], L0=config['L'], H=config['H'], Nh=2**config['LEVEL'], Nl=config['NL'])
    # Assemble into xarray dataset 
    ds = xr.Dataset(data_vars={name: (['t','zl','x','y'], np.expand_dims(array, axis=0)) for name, array in zip(fieldnames, fields)},
                    coords=dict(t=(['t'], np.array([t])),
                                x=(['x'], 0.5*(x_mesh[0,:-1,0]+x_mesh[0,1:,0])), 
                                x_g=(['x_g'], x_mesh[0,:-1,0]),
                                y=(['y'], 0.5*(y_mesh[0,0,:-1]+y_mesh[0,0,1:])),
                                y_g=(['y_g'], y_mesh[0,0,:-1]),
                                zl=(['zl'], np.arange(0,config['NL'])),
                                zl_g=(['zl_g'], np.arange(0,config['NL'])-0.5),
                                z=(['t','zl','x','y'], 
                                   np.expand_dims(0.5*(z_mesh[:-1,:-1,:-1]+z_mesh[1:,:-1,:-1]), axis=0)),
                                z_g=(['t','zl_g','x','y'], np.expand_dims(z_mesh[:-1,:-1,:-1], axis=0))),
                    attrs=dict(sourcepath=field_path, **config))
    # Reading eta (seems not necessary) since eta can be retrieved from sum of h
    # And there is a slight time shift if vorticity file is generated later than eta
#     matrix = np.fromfile(eta_path+'eta_matrix_%g' %t, dtype=np.float32)
#     N = 2**config['LEVEL']; matrix = matrix.reshape(N+1,N+1); matrix = matrix[1:,1:]
#     ds['eta'] = (['x','y'], matrix)

    # Same as float16 to save space? Optional. Can save about half the space.
    ds = ds.astype('float32')
    return ds

''' Pickle objects. '''
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)
def load_object(filename):
    with open(filename, 'rb') as input:  # Overwrites any existing file.
        obj = pickle.load(input)
    return obj