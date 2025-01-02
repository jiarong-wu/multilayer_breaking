import numpy as np

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

''' An old way of different different time from different folders '''
def read_t_old(self, t=5):
    h_ensem = []; ux_ensem = []; uy_ensem = []; uz_ensem = []
    folder = self.path + 'field_t%g/' %t
    N = 2**self.LEVEL
    for l in range (0,self.NL):
        h = read(filepath=folder, pre='h_matrix_l', index=l, N=N)
        ux = read(filepath=folder, pre='ux_matrix_l', index=l, N=N)
        uy = read(filepath=folder, pre='uy_matrix_l', index=l, N=N)
        uz = read(filepath=folder, pre='uz_matrix_l', index=l, N=N)
        h_ensem.append(h)
        ux_ensem.append(ux)
        uy_ensem.append(uy)
        uz_ensem.append(uz)
    """ axis0-z; axis1-x; aixs2-y"""
    self.h_ensem = np.array (h_ensem)
    self.ux_ensem = np.array (ux_ensem)
    self.uy_ensem = np.array (uy_ensem)
    self.uz_ensem = np.array (uz_ensem)
    return (self.h_ensem, self.ux_ensem, self.uy_ensem, self.uz_ensem)