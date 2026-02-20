import numpy as np
import warnings
warnings.filterwarnings("ignore")

from scipy.stats import rv_continuous

class trunc_isotropic_sin(rv_continuous):
    """Truncated isotropic(sin) distribution
    """
    name = "trunc_isotropic_sin"
    def __init__(self,x_low,x_high):
        
        self.x_low = x_low
        self.x_high = x_high
        super(trunc_isotropic_sin, self).__init__(a=x_low, b=x_high)

    def _pdf(self,x):
        """ Probability density function of the truncated isotropic sin distribution
        """
        		
        k = np.sin(x)/(np.cos(self.x_low) - np.cos(self.x_high))

        return k