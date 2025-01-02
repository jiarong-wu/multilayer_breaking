
""" The detection function """
import os
import numpy as np
from skimage.feature import hessian_matrix, hessian_matrix_eigvals

def detect_ridges(gray, sigma=1):
    H_elems = hessian_matrix(gray, sigma=sigma, order='rc')
    maxima_ridges, minima_ridges = hessian_matrix_eigvals(H_elems)
    return maxima_ridges, minima_ridges

def detect_slope(grad, threshold):
    return (np.logical_not(grad > threshold))
    
    
""" 
Define the class Case (for each realization, varying N spreading, Nx, Ny, Nl etc.) and the overarching class Config (for the same L0 and P).
(What about separating the definition of Cases).
Cases contains ridge detection and simple mapping methods. """

class Case(object):
    """ This class defines methods specific to cases.
        Attributes: 
            self.path
            self.config (L0,P)
            self.level        
        
            self.hist_ensem, self.bins (by calling simple_mapping)
    """
    
    def read_t(self, fieldname='eta', t=0):
        filename = self.path + 'surface/'+ fieldname + '_matrix_%g' %t
        isFile = os.path.isfile(filename)
        if isFile != True:
            print (filename + ' does not exist!') 
            return 1
        matrix = np.fromfile(filename, dtype=np.float32)
        N = 2**self.LEVEL
        matrix = matrix.reshape(N+1,N+1); matrix = matrix[1:,1:]
        matrix = np.rot90(matrix)
        return matrix
    
    def simple_mapping(self, tseries, bins=[], method=0, threshold=0, EXTRA_FILTER=True):
        """ Arguments: 
                tseries: the time where the stats are collected 
                bins: velocity bins where the crest is counted
                method: 0-ridge detection; 1-slope thresholding
                threshold: if method=0, an empirical (ad hoc) threshold to filter the ridges; if method=1, the slope threshold
            Outputs: 
                hist_taver: 1D array of velocity binned crest length (in unit of pixels), averaged over tseries.  
        """    
        """ Iteratively read fields and extract crests """ 
        hist_t = []
        N = 2**self.LEVEL
        for t in tseries:
            eta = self.read_t('eta', t)
            ux = self.read_t('ux', t)
            uy = self.read_t('uy', t)

            """ Edge detection """
            if method == 0:
                a, b = detect_ridges(eta, sigma=1.0) # Maxima and minima ridges
                delta = self.config.L0/2**self.LEVEL # Normalize the curvature by grid size
                b_norm = b/delta**2                
                b_ = np.logical_not(b_norm > threshold) # Is this value fixed???
                # Extra filtering by above 2.5 sigma (we can kind of use something between 2.5-3\sigma)
                if EXTRA_FILTER == True:
                    height_filter = np.logical_not(eta < 2.5*np.var(eta)**0.5)
                    b_ = b_*height_filter
            elif method == 1:
                b_ = detect_slope(eta_gradx, threshold=threshold)         
            """ One sided edge """
            a_ = np.zeros(b_.shape) 
            for i in range(0,N-1):
                for j in range(1,N-1):
                    if (b_[i][j-1] > 0) and (b_[i][j+1] == 0):
                        a_[i][j] = 1
                        
            # Extra filtering by above 2.5 sigma (we can kind of use something between 2.5-3\sigma)
#             if EXTRA_FILTER == True:
#                 height_filter = np.logical_not(eta < 2.5*np.var(eta)**0.5)
#                 a_ = a_*height_filter 
                
            ux_a = ux*a_
            uy_a = uy*a_ 
            hist,_ = np.histogram((ux_a**2+uy_a**2)**0.5, bins=bins, density=False)
            hist_t.append(hist)   
            hist_taver = np.average(np.array(hist_t), axis=0)
        return hist_taver
    
    def time_window(self, tseries_ensem, bins=[], method=0, threshold=0):
        """ This function collects stats with a series of time windows
            Arguments: 
                tseries_ensem: an list of time windows 
            Outputs:
                self.hist_ensem: a dictionary with tseries_ensem and the corresponding hist_ensem (with shared crest velocity bins and threshold)
                                 self.hist_ensem = {"tseries_ensem": tseries_ensem, "hist_ensem": hist_ensem, "bins", "threshold": threshold}
        """
        """ Box size dependent threshold and bins """
        if threshold == 0: # if not specified: use the default one based on box size
#             if self.config.L0 == 50:
#                 threshold = -0.006
#             elif self.config.L0 ==200:
#                 threshold = -0.1
#             elif self.config.L0 == 500:
#                 threshold = -0.07
#             else:
#                 print('Need to specify threshold!')
#                 return 1
            threshold = -3*self.config.kp

        if bins == []:
            if self.config.kp == 2*np.pi/10:
                bins = np.array([0.01])
                bins = np.concatenate((bins,np.arange(0.25,5,0.25)))   
            elif self.config.kp == 2*np.pi/40:
                bins = np.array([0.01])
                bins = np.concatenate((bins,np.arange(0.5,10,0.25)))   
            elif self.config.kp == 2*np.pi/100:
                bins = np.array([0.01])
                bins = np.concatenate((bins,np.arange(0.5,15,0.25)))   
            else:
                print('Need to specify bins!')
                return 1 
        hist_ensem = []           
        for t_series in tseries_ensem:
            hist_taver = self.simple_mapping(t_series, bins, method, threshold)        
            hist_ensem.append(hist_taver)
        self.hist_ensem = {"tseries_ensem": tseries_ensem, "hist_ensem": hist_ensem, "bins":bins, "threshold": threshold}
        


""" Define the mega-data of cases """
class Config():
    """ Input:
            L0, P, kp: A few parameters
            Prefix: either 'field_new_' or' field_new_500m_'
    """

    def __init__(self, L0=50, P=0.01, kp=0, prefix='field_new_'):
        self.cases = [] # Empty list of cases  
        self.L0 = L0; self.P = P; self.kp = kp; self.prefix=prefix
        if self.kp == 0:
            self.kp = 2*np.pi/self.L0*5
        self.g = 9.8
        self.cp = (self.g/self.kp)**0.5      
                  
    def add_realization(self, LEVEL=9, NL=30, rand=0, RE=40000, Htheta=0.5, Npower=5, prefix='', postfix='noHtheta/', path=''):
        """ Call this method to define realizations of the same configuration. In this way all the cases that belong to the same configuration
            are grouped together. 
            The Realization class inherites the Case class (defined outside Config) but in addition it knows about the overarching configuration. """
        config = self # this is the outer self
        class Realization(Case):
            def __init__(self): # this is the inner self
                self.config = config
                self.LEVEL = LEVEL; self.NL = NL; self.rand = rand; self.RE = RE; self.Htheta=Htheta; self.Npower=Npower
                self.path = self.config.prefix + prefix + 'P%g' %self.config.P + '_RE%g_' %RE + '%g_' %LEVEL + '%g_' %NL + 'rand%g_' %rand + postfix
                if path != '':                    
                    if path != self.path:
                        print('Given path does not agree with parameter values!') # warning but no error
                        self.path = path
                isFile = os.path.isfile(self.path)
                if isFile != 0:
                    print ('Folder does not exist!') 
                    return 1
                config.cases.append(self)
        return Realization()
    
    def delete_realization(self, index=-1):
        self.cases.remove(self.cases[index])
        print('%g th case removed, %g left' %(index, len(self.cases)))
    