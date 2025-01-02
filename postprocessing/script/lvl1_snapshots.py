import argparse
import re
import os 
from mlpython.fio import read_netcdf
import gc

##### Change these if needed before running the python file #####
# tseries: the times 
# fieldnames: variable names to read
# savepath: the parent path for all processed cases
# Example:
# python read.py --path='/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.02_RE40000_10_15_rand2_Htheta0.503/' --label='C4'
tseries = (0,40,80,120,160,110,130,150,170)
fieldnames = ['h','ux','uy','uz']
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
dirpath = savepath + config['label']
os.makedirs(dirpath, exist_ok=True)  # `exist_ok=True` prevents error if directory exists
print(f"Directory '{dirpath}' created!")

for t in tseries:
    print('Reading t=%g...' %t)
    ds = read_netcdf (args.path + 'vort/', config, t, fieldnames)
    filename = dirpath + '/field%g.nc' %t  
    encoding = {}
    for var_name in ds.data_vars:
        encoding[var_name] = {'dtype': 'float32', 'zlib': True}
    ds.to_netcdf(filename, encoding=encoding)
    
    del ds # this is needed for memory 
    gc.collect()