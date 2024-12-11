import xarray as xr

''' Function to read files and assemble them into an xarray dataset '''
def read_netcdf (path, config, t, fieldnames=['h','ux','uy','uz','omegax','omegay','omegaz','dzdx','dzdy','dzdxc','dzdyc']):
    # 3D field path and 2D eta path relative to the main path
    field_path = path + 'vort/'
    eta_path = path + 'surface/'
    # Reading fields (centered) according to the name list
    fields = read_t(fieldnames=fieldnames, t=t, Nh=2**config['LEVEL'], Nl=config['NL'], path=field_path)    
    # 3D meshes at vertices
    x_mesh, y_mesh, z_mesh = array_to_mesh (fields[fieldnames.index('h')], L0=config['L'], H=config['H'], Nh=2**config['LEVEL'], Nl=config['NL'])
    # Assemble into xarray dataset 
    ds = xr.Dataset(data_vars={name: (['t','zl','x','y'], np.expand_dims(array, axis=0)) for name, array in zip(fieldnames, fields)},
                    coords=dict(x=(['x'], 0.5*(x_mesh[0,:-1,0]+x_mesh[0,1:,0])), 
                                x_g=(['x'], x_mesh[0,:-1,0]),
                                y=(['y'], 0.5*(y_mesh[0,0,:-1]+y_mesh[0,0,1:])),
                                y_g=(['y'], y_mesh[0,0,:-1]),
                                zl=(['zl'], np.arange(0,config['NL'])),
                                zl_g=(['zl'], np.arange(0,config['NL'])-0.5),
                                z=(['t','zl','x','y'], 
                                   np.expand_dims(0.5*(z_mesh[:-1,:-1,:-1]+z_mesh[1:,:-1,:-1]), axis=0)),
                                z_g=(['t','zl','x','y'], np.expand_dims(z_mesh[:-1,:-1,:-1], axis=0))),
                    attrs=dict(sourcepath=path, **config))
    # Reading eta (seems not necessary) since eta can be retrieved from sum of h
    # And there is a slight time shift if vorticity file is generated later than eta
#     matrix = np.fromfile(eta_path+'eta_matrix_%g' %t, dtype=np.float32)
#     N = 2**config['LEVEL']; matrix = matrix.reshape(N+1,N+1); matrix = matrix[1:,1:]
#     ds['eta'] = (['x','y'], matrix)

    # Same as float16 to save space? Optional. Can save about half the space.
    ds = ds.astype('float32')
    return ds