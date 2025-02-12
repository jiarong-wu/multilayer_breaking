''' Compute the breaking stats ''' 
import numpy as np
import xarray as xr
import gc
from mlpython.breaking import simple_mapping, get_bins
from dask.diagnostics import ProgressBar

# paths = ['/projects/DEIKE/jiarongw/multilayer/JPO/processed/C1',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C2',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C3',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C4',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C5',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C4',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C4_rand4',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C5',
#          '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C5_rand4',]
#### Files that admits finer time stepping, for smoother stats lines ####
paths = ['/projects/DEIKE/jiarongw/multilayer/JPO/processed/C1',
         '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C2',
         '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C3',
         '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C4_rand4',
         '/projects/DEIKE/jiarongw/multilayer/JPO/processed/C5']

paths = ['/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_NL30',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_NL45',]

#### Pick time windows of analyzing ####
time = np.arange(100,181,1)

#### Parallel computing doesn't seem to work yet ####
# from dask.distributed import Client
# client = Client(processes=False)  # threads only
# print(client.dashboard_link)  # View the dashboard at the provided link

for path in paths:
    ds = xr.open_dataset(path + '/series.nc', chunks={'t': 100, 'x':-1, 'y':-1})
    ds = ds.sel(t=time, method='nearest')
    #### Some metadata ####
    delta = ds.attrs['L']/2**ds.attrs['LEVEL'] # Normalize the curvature by grid size
    kp = 2*np.pi/ds.attrs['L']*5
    threshold = -3*kp
    bins = get_bins(kp)
    bins_center = bins[1:] - (bins[2] - bins[1])/2
    #### Compute stats ####
    print('Computing stats for %s...' %path)
    with ProgressBar():
        hist = xr.apply_ufunc (
            simple_mapping,
            ds.eta, 
            ds.ux,  
            ds.uy,
            input_core_dims=[['x','y'], ['x','y'], ['x','y']],  # Core dimensions for each input
            output_core_dims=[['c']],  # Core dimensions for the output
            # exclude_dims=set(('zl',)),
            output_sizes={'c':len(bins_center)},
            output_dtypes=['float32'],
            vectorize=True,  # Enable vectorization
            dask="parallelized",  # Parallelize using Dask if the data is large
            kwargs={'delta':delta, 'bins':bins, 'method':0, 'threshold':threshold} 
        ).compute()
    hist = hist.assign_coords(c=bins_center)
    
    #### Save to a separate file ####
    filename = path + '/breaking_hist.nc'
    compression_settings = {
        'zlib': True,
        'complevel': 5  # Compression level from 1 (fastest, least compression) to 9 (slowest, most compression)
    }
    hist.name = 'hist'
    hist.to_netcdf(filename, encoding={'hist': compression_settings})
    print('Breaking stats saved!')
    del(hist) # Delete ds for memory 
    gc.collect()

# """ Compute the energy loss (without filtering) and the breaking stats """
# for k, config in enumerate(config_set[:5]):
#     for case in config.cases:
#         if (case.NL == 15) and (case.LEVEL == 10) and (case.rand != 0) and (case.rand != 1) and (case.Npower == 5):
#             print (case.path)
#             energy = pd.read_table(case.path +'energy_after_remap.dat', delimiter=' ', names=['t','ke','gpe'])
#             energy = energy.drop_duplicates(subset=['t'])
#             case.dEkdt = []; case.dEpdt = [] 
#             for tseries in tseries_ensem:
#                 print('From t = %g to %g.' %(tseries[0], tseries[-1]))
#                 idx1 = (np.abs(energy.t - tseries[0])).argmin()
#                 idx2 = (np.abs(energy.t - tseries[-1])).argmin()
#                 dEkdt = (energy['ke'].values[idx1] - energy['ke'].values[idx2])/(tseries[-1]-tseries[0])
#                 dEpdt = (energy['gpe'].values[idx1] - energy['gpe'].values[idx2])/(tseries[-1]-tseries[0])           
#                 case.dEkdt.append(dEkdt); case.dEpdt.append(dEpdt)
#             case.dEkdt = np.array(case.dEkdt); case.dEpdt = np.array(case.dEpdt)
#             case.dEdt = (case.dEkdt + case.dEpdt)/case.config.L0**2 # dissipation per area
#             print(case.dEdt)
             
#             picklename = case.path + 'breakingstat.pkl'   
#             ''' If restore from pickle '''
#             if os.path.isfile(picklename):
#                 case.hist_ensem = load_object(picklename)
#             else:
#                 case.time_window(tseries_ensem, threshold=0, bins=[])
#                 save_object(case.hist_ensem, picklename)
                
# """ Pickle the breaking statistics """
# for k, config in enumerate(config_set[:4]):
#     for case in config.cases:
#         if (case.NL == 15) and (case.LEVEL == 10) and (case.rand != 0) and (case.rand != 1) and (case.Npower == 5):
#             picklename = case.path + 'breakingstat.pkl'
#             save_object(case.hist_ensem, picklename)