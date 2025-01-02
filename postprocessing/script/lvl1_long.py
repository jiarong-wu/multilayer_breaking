import argparse
import re
from mlpython.fio import read_netcdf

##### Change these if needed before running the python file #####
# tseries: the times 
# fieldnames: variable names to read
# savepath: the parent path for all processed cases
tseries = (0,40,80,120,160,110,130,150,170)
fieldnames = ['h','ux','uy','uz','omegax','omegay','omegaz']
savepath = '/projects/DEIKE/jiarongw/multilayer/JPO/processed/'

if __name__ == "__main__":
    doc = 'Read the raw output from runs.'
    parser = argparse.ArgumentParser(description=doc)
    parser.add_argument('--path', type=str, required=True, 
                        help='Folder path (absolute).')
    parser.add_argument('--label', type=str, required=True, 
                        help='Case label (used to name processed folders).')
    parser.add_argument('--t', type=float, required=True,
                        help='Simulation time.')
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
    ds = read_netcdf (args.path, config, fieldnames)
    filename = savepath + config['label'] + '/field%g.nc' %t
    encoding = {}
    for var_name in ds.data_vars:
        encoding[var_name] = {'dtype': 'float32', 'zlib': True}
    ds.to_netcdf(filename, encoding=encoding)
    
#### Pick time windows of analyzing ####
tseries_ensem = []
for i in range(0,4):
    tstart = 100 + 20*i
    dt = 0.2
    tseries = np.arange(tstart, tstart+20, dt)
    tseries_ensem.append(tseries)
    

""" Compute the energy loss (without filtering) and the breaking stats """
for k, config in enumerate(config_set[:5]):
    for case in config.cases:
        if (case.NL == 15) and (case.LEVEL == 10) and (case.rand != 0) and (case.rand != 1) and (case.Npower == 5):
            print (case.path)
            energy = pd.read_table(case.path +'energy_after_remap.dat', delimiter=' ', names=['t','ke','gpe'])
            energy = energy.drop_duplicates(subset=['t'])
            case.dEkdt = []; case.dEpdt = [] 
            for tseries in tseries_ensem:
                print('From t = %g to %g.' %(tseries[0], tseries[-1]))
                idx1 = (np.abs(energy.t - tseries[0])).argmin()
                idx2 = (np.abs(energy.t - tseries[-1])).argmin()
                dEkdt = (energy['ke'].values[idx1] - energy['ke'].values[idx2])/(tseries[-1]-tseries[0])
                dEpdt = (energy['gpe'].values[idx1] - energy['gpe'].values[idx2])/(tseries[-1]-tseries[0])           
                case.dEkdt.append(dEkdt); case.dEpdt.append(dEpdt)
            case.dEkdt = np.array(case.dEkdt); case.dEpdt = np.array(case.dEpdt)
            case.dEdt = (case.dEkdt + case.dEpdt)/case.config.L0**2 # dissipation per area
            print(case.dEdt)
             
            picklename = case.path + 'breakingstat.pkl'   
            ''' If restore from pickle '''
            if os.path.isfile(picklename):
                case.hist_ensem = load_object(picklename)
            else:
                case.time_window(tseries_ensem, threshold=0, bins=[])
                save_object(case.hist_ensem, picklename)
            
# """ Pickle the breaking statistics """
# for k, config in enumerate(config_set[:4]):
#     for case in config.cases:
#         if (case.NL == 15) and (case.LEVEL == 10) and (case.rand != 0) and (case.rand != 1) and (case.Npower == 5):
#             picklename = case.path + 'breakingstat.pkl'
#             save_object(case.hist_ensem, picklename)
    
    