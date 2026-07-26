import numpy as np
import warnings
warnings.filterwarnings("ignore")

from scipy.stats import rv_continuous
from scipy.special import gamma,gammaincc


class trunc_gamma(rv_continuous):
    """Truncated gamma distribution
    """
    name = "trunc_gamma"
    def __init__(self,scale,shape,x_low,x_high):
        
        self.scale = scale
        self.shape = shape
        self.x_low = x_low
        self.x_high = x_high
        super(trunc_gamma, self).__init__(a=x_low, b=x_high)

    def _pdf(self,x):
        """ Probability density function of the truncated gamma distribution
        
        see https://arxiv.org/pdf/1401.0287.pdf
        see https://arxiv.org/pdf/1912.12053.pdf
        """
        
        den = self.scale*gamma(1+self.shape)*(gammaincc(1+self.shape,self.x_low/self.scale) -gammaincc(1+self.shape,self.x_high/self.scale)) \
            +np.exp(-self.x_high/self.scale)*self.scale**(-self.shape+1)*self.x_high**self.shape - np.exp(-self.x_low/self.scale)*self.scale**(-self.shape+1)*self.x_low**self.shape
        num = self.shape 
        k = num/den
        return k*(x/self.scale)**(self.shape-1)*np.exp(-x/self.scale)