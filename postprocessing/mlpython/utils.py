from scipy.interpolate import interp1d
import xarray as xr

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