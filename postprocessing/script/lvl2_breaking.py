''' Compute the breaking stats ''' 

paths = ['/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_rand4',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C1',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C5_rand4',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C3',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C5',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C2',
         '/Users/jiarongw/Data/multilayer_data/JPO2024/processed/C4_NL45']


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