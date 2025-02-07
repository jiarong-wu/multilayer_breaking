''' 
This script compute horizontally averaged dissipation across different cases 
and save them to one netcdf file. 
'''
import os
import gc
import xarray as xr
import numpy as np
from mlpython.utils import compute_diss
from mlpython.utils import interpz


# Suppress all warnings
import warnings
warnings.filterwarnings("ignore")

#### Data path ####
base_dir = "/Users/jiarongw/Data/multilayer_data/JPO2024/processed/"
#### List of time ####
tlist = np.array([110,130,150,170])
#### Interpolation grid ####
znew = np.arange(-20,1,0.1)

#### Find all cases ####
# paths = []
# for item in os.listdir(base_dir):
#     full_path = os.path.join(base_dir, item)
#     paths.append(full_path)
paths = ['/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_rand4',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C1',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C5_rand4',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C3',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C5',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C2',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_NL30',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_NL45']

paths = ['/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_NL30',]

for path in paths:
    print('Reading... dir='+path)
    filelist = [path + f'/field{t}.nc' for t in tlist]
    ds = xr.open_mfdataset(filelist, concat_dim='t', combine='nested')
    ds = ds.assign_coords(t=tlist)
 
    # Compute dissipation
    ds = compute_diss(ds).compute()
    
    # Save 1D epsilon
    ds1d = xr.Dataset(data_vars={'epsilon':(['t','zl'], ds.epsilon.mean(dim=['x','y']).values),
                                'omegax':(['t','zl'], ds.epsilon.mean(dim=['x','y']).values),
                                'omegay':(['t','zl'], ds.epsilon.mean(dim=['x','y']).values),
                                'omegaz':(['t','zl'], ds.epsilon.mean(dim=['x','y']).values),},
                    coords={'t':(['t'], ds.t.values),
                            'z':(['t','zl'], ds.z.mean(dim=['x','y']).values)},
                    attrs={'description':'Layer averaged epsilon and vorticity.', **ds.attrs})

    filename = path + '/epsilon1d.nc'
    encoding = {}
    for var_name in ds1d.data_vars:
        encoding[var_name] = {'dtype': 'float32', 'zlib': True}
    ds1d.to_netcdf(filename, encoding=encoding, engine='h5netcdf')
    print('Layer averaged 1D profiles saved!')

    # Interpolate to cartesian grid
    fields = ['epsilon','omegaxp','omegayp','omegazp']
    fields_t = []
    for field in fields:
        fieldt = []
        for t in ds.t:
            field_interp = interpz(znew, ds.z.sel(t=t), ds[field].sel(t=t), fill_value=np.nan).mean(dim=['x','y']).compute() 
            fieldt.append(field_interp)
        fields_t.append(np.array(fieldt))
        
    ds1d_interp = xr.Dataset(data_vars={name: (['t','z'], array) for name, array in zip(fields, fields_t)},
                            coords={'t':(['t'], ds.t.values),
                                    'z':(['z'], znew)},
                            attrs={'description':'Interpolated and averaged epsilon and vorticity.', **ds.attrs})

    filename = path + '/epsilon1d_interp.nc'
    encoding = {}
    for var_name in ds1d_interp.data_vars:
        encoding[var_name] = {'dtype': 'float32', 'zlib': True}
    ds1d_interp.to_netcdf(filename, encoding=encoding, engine='h5netcdf')
    print('Interpolated averaged 1D profiles saved!')
    
    # Delete ds for memory 
    del(ds, ds1d, ds1d_interp)
    gc.collect()

