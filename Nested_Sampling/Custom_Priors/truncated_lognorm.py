import numpy as np
import warnings
warnings.filterwarnings("ignore")

from scipy.stats import rv_continuous
from scipy.special import erf


class trunc_lognorm(rv_continuous):
    """Truncated lognormal distribution
    """
    name = "trunc_lognorm"
    def __init__(self,scale,shape,x_low,x_high):
        
        self.scale = scale
        self.shape = shape
        self.x_low = x_low
        self.x_high = x_high
        super(trunc_lognorm, self).__init__(a=x_low, b=x_high)

    def _pdf(self,x):
        """ Probability density function of the truncated lognormal distribution
        
        see https://arxiv.org/pdf/1708.06159.pdf
        """
        
        den = -np.sqrt(np.pi)*self.shape*(erf(0.5*(np.sqrt(2)/self.shape)*np.log(self.x_low/self.scale)) - erf(0.5*(np.sqrt(2)/self.shape)*np.log(self.x_high/self.scale)))*(x)
        num = np.sqrt(2)*np.exp(-0.5*(1/self.shape**2)*(np.log((x)/self.scale))**2)

        k = num/den
        k = np.nan_to_num(k)

        return k