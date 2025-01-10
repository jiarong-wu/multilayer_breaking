# ''' Compute the breaking stats ''' 
# import numpy as np
# import xarray as xr
# import dask.array as da
# import gc
# from dask.diagnostics import ProgressBar

# # Create a synthetic DataArray
# data = np.random.rand(10000, 10000)
# hist = xr.DataArray(data, dims=['x', 'y'])

# # Convert the DataArray to a Dask array
# hist = hist.chunk({'x': 100, 'y': 100})  # Adjust chunk size as needed

# # Define a simple custom function
# def custom_function(x):
#     return x ** 2

# # Apply the custom function in parallel using Dask
# result = xr.apply_ufunc(
#     custom_function,
#     hist,
#     dask='parallelized',
#     output_dtypes=[hist.dtype]
# )

# # Save the result to a NetCDF file with compression
# filename = '/projects/DEIKE/jiarongw/multilayer/JPO/processed/breaking_hist.nc'
# compression_settings = {
#     'zlib': True,
#     'complevel': 5  # Compression level from 1 (fastest, least compression) to 9 (slowest, most compression)
# }
# result.name = 'hist'

# # Use a progress bar to monitor the computation
# with ProgressBar():
#     result.to_netcdf(filename, encoding={'hist': compression_settings})

# print('Breaking stats saved!')
# del(result)  # Delete ds for memory
# gc.collect()

import torch
import torch.multiprocessing as mp
import gc

def simple_operation(tensor):
    return tensor ** 2

def worker(rank, tensor, results, device):
    tensor = tensor.to(device)
    result = simple_operation(tensor)
    results[rank] = result.cpu()
    print(f'Worker {rank} finished on device {device}')

if __name__ == '__main__':
    # Check if CUDA is available and set the device
    if torch.cuda.is_available():
        devices = [torch.device(f'cuda:{i}') for i in range(torch.cuda.device_count())]
    else:
        devices = [torch.device('cpu')]

    # Create a synthetic tensor
    tensor = torch.rand(100000, 100000)

    # Split the tensor for parallel processing
    split_tensors = torch.chunk(tensor, len(devices))

    # Create a shared list to store results
    manager = mp.Manager()
    results = manager.list([None] * len(devices))

    # Start the parallel processing
    processes = []
    for rank, device in enumerate(devices):
        p = mp.Process(target=worker, args=(rank, split_tensors[rank], results, device))
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    # Combine the results
    result = torch.cat(tuple(results))

    # Save the result to a file
    # torch.save(result, '/projects/DEIKE/jiarongw/multilayer/JPO/processed/result.pt')
    # print('Result saved!')

    del(result)
    gc.collect()