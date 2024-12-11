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
    
    