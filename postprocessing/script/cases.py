paths = ['/projects/DEIKE/jiarongw/multilayer/JFM/field_new_200m_P0.008_RE40000_10_15_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/JFM/field_new_200m_P0.01_RE40000_10_15_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/JFM/field_new_200m_P0.016_RE40000_10_15_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.02_RE40000_10_15_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.02_RE40000_10_15_rand4_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.02_RE40000_10_30_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.02_RE40000_10_45_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.03_RE40000_10_15_rand2_Htheta0.503/',
         '/projects/DEIKE/jiarongw/multilayer/revision/field_new_200m_P0.03_RE40000_10_15_rand4_Htheta0.503/',]

labels = ['C1','C2','C3','C4','C4_rand4','C4_NL30','C4_NL45','C5','C5_rand4']

savepath = '/projects/DEIKE/jiarongw/multilayer/JPO/processed/'

import os
for path, label in zip(paths, labels):
    os.system(f'cp {path}/energy_before_remap.dat {savepath}{label}/energy.dat')