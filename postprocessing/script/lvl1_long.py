''' surface/ folder has eta, ux, uy fields with temporal resolution of 0.1.
    energy_before_remap.dat file has integrated energy with non-uniform temporal resolution (used i++) around 0.02.
'''
import argparse
import re
import numpy as np
import xarray as xr
import pandas as pd

def read_long (field_path, config, time, fieldnames=['eta','ux','uy']):
    #### Define time window, sampling interval, and a few metadata ####
    # field_path = '/projects/DEIKE/jiarongw/multilayer/JFM/field_new_200m_P0.008_RE40000_10_15_rand2_Htheta0.503/' 
    Nt = len(time)
    N = 2**config['LEVEL'] # Change to config  
    L0 = config['L'] # Change to config

    #### Read in 2D fields ####
    fieldnames = ['eta','ux','uy']
    fields = []
    for fieldname in fieldnames:
        f_series = np.zeros((Nt,N,N), dtype=np.float32)
        for i in range(0, Nt):
            filename = field_path + f'surface/{fieldname}_matrix_%g' %time[i]
            f = np.fromfile(filename, dtype=np.float32)
            f = f.reshape(N+1,N+1); f = f[1:,1:]
            f_series[i] = f  
        fields.append(f_series)    
        
    x_mesh = np.linspace(-L0/2, L0/2, N+1, endpoint=True)
    ds = xr.Dataset(data_vars={name: (['t','x','y'], array) for name, array in zip(fieldnames, fields)},
                    coords={'t': (['t'], time),
                            'x': (['x'], 0.5*(x_mesh[:-1]+x_mesh[1:])),
                            'y': (['y'], 0.5*(x_mesh[:-1]+x_mesh[1:]))},
                    attrs=dict(sourcepath=field_path, **config))

    #### Read in energy and interpolate onto surface field time ####
    energy = pd.read_table(field_path +'energy_before_remap.dat', delimiter=' ', names=['t','ke','gpe'])
    energy = energy.drop_duplicates(subset=['t'])
    df = energy.set_index('t'); ds_energy = df.to_xarray()
    ds['ke'] = ds_energy.ke.interp(t=ds.t)
    ds['gpe'] = ds_energy.gpe.interp(t=ds.t)
    
    return ds

##### Change these if needed before running the python file #####
# time: a long series with small interval
# fieldnames: variable names to read (surface fields)
# savepath: the parent path for all processed cases
time = np.arange(100,181,1) # For some cases eta saved every t=1
# time = np.arange(100,120.1,0.1) # For others 0.1
fieldnames = ['eta']
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
    print('Assembling...')
    ds = read_long (args.path, config, time, fieldnames)
    filename = savepath + config['label'] + '/series.nc' 
    encoding = {}
    for var_name in ds.data_vars:
        encoding[var_name] = {'dtype': 'float32', 'zlib': True}
    ds.to_netcdf(filename, encoding=encoding)
    print('Long series of surface fields and energy saved to ' + filename)