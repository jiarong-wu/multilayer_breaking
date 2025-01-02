import xarray as xr
import argparse
import re
from bulk_func import array_to_mesh, read_t
import numpy as np
import os
# from mlpython.fio import read_netcdf

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
                    coords=dict(x=(['x'], 0.5*(x_mesh[0,:-1,0]+x_mesh[0,1:,0])), 
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


##### Change these if needed before running the python file #####
# tseries: the times 
# fieldnames: variable names to read
# savepath: the parent path for all processed cases
# Example:
# python read.py --path='/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.02_RE40000_10_15_rand2_Htheta0.503/' --label='C4'
tseries = (80,0)
fieldnames = ['h','ux','uy','uz','omegax','omegay','omegaz']
savepath = '/projects/DEIKE/jiarongw/multilayer/JPO/processed/'

if __name__ == "__main__":
    doc = 'Read the raw output from runs.'
    parser = argparse.ArgumentParser(description=doc)
    parser.add_argument('--path', type=str, required=True, 
                        help='Folder path (absolute).')
    parser.add_argument('--label', type=str, required=True, 
                        help='Case label (used to name processed folders).')
    args = parser.parse_args()
    
# Extract parameters from path string
numbers = re.findall(r"-?\d+\.?\d*", args.path)
numbers = [float(num) if '.' in num else int(num) for num in numbers]

# Assemble config dict
names = ['L', 'P', 'Re', 'LEVEL', 'NL', 'rand', 'Htheta']
config = {name: number for name, number in zip(names, numbers)}
config = {**config, 'H': config['L']/5, 'label':args.label} # Add the depth that is missing in the folder name

# Read and save to individual netcdf for each time
for t in tseries:
    print('Reading t=%g...' %t)
    ds = read_netcdf (args.path + 'vort_new_new_new/', config, t, fieldnames)
    
    dirpath = savepath + config['label']
    os.makedirs(dirpath, exist_ok=True)  # `exist_ok=True` prevents error if directory exists
    print(f"Directory '{dirpath}' created!")
    filename = dirpath + '/field%g.nc' %t
    
    encoding = {}
    for var_name in ds.data_vars:
        encoding[var_name] = {'dtype': 'float32', 'zlib': True}
    ds.to_netcdf(filename, encoding=encoding)
    