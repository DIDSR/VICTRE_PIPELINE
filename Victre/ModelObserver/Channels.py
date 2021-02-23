from .LinSpace import LinSpace
from .Exceptions import checkChannelName, checkChannelArgs, defaultParams
from collections import deque
from scipy.special import comb, factorial
import numpy as np
import math


class Channels:
    def __init__(self, dimensions, channel_params, ftype):
        '''
        Initialize the 2D channels for the CHO model.

        Parameters
        ----------
        dimensions : TYPE, tuple of length 2.
            DESCRIPTION. Dimensions of the image (rows,columns)

        Returns
        -------
        None.

        '''
        self.dim = dimensions
        self.ftype = ftype
        self.channel_params = channel_params

    def _norm_channels(self, channel):
        '''
        Normalize the values in the 2D channel template. 

        Parameters
        ----------
        channel : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.

        '''
        return (channel / np.sum(channel))

    def _frequency_DoG(self, n, q=1.67, alpha=1.4, s0=.005, norm=True):
        '''
        Create DoG filter in the frequency domain??

        Parameters
        ----------
        n : TYPE, integer.
            DESCRIPTION. 
        q : TYPE, optional
            DESCRIPTION. The default is 1.67.
        alpha : TYPE, optional
            DESCRIPTION. The default is 1.4.
        s0 : TYPE, optional
            DESCRIPTION. The default is .005.

        Returns
        -------
        DoG filter in the frequency domain with the DC component in the 
        center and a dictionary of the parameters used to create the DoG.

        '''
        x, y = self._euclidean2D()
        u, v = x / self.dim[1], y / self.dim[0]
        freqs = np.sqrt(u**2 + v**2)  # 2D spatial frequency domain

        sj = s0 * (alpha**n)
        g1 = np.exp((-1 / 2) * (freqs / (q * sj))**2)
        g2 = np.exp((-1 / 2) * (freqs / sj)**2)

        # params = {'n': n,
        #           'q': q,
        #           'alpha': alpha,
        #           's0': s0,
        #           }
        dog = g1 - g2
        if norm is True:
            dog = self._norm_channels(dog)
        return dog

    def _spatial_Laguerre_Gaussian(self, a=30, b=25, n=5, norm=True):
        '''
        Generate a Laguerre—Gauss Filter in the spatial domain. 

        Parameters
        ----------
        a : TYPE, optional
            DESCRIPTION. The default is 8.
        b : TYPE, optional
            DESCRIPTION. The default is 8.
        n : TYPE, optional
            DESCRIPTION. The default is 6.

        Returns
        -------
        TYPE, 2D numpy array
            DESCRIPTION. The filter in the spatial domain.  

        '''
        xc, yc = 0, 0
        x, y = self._euclidean2D()

        def lgpoly(x):
            total = 0
            for m in range(n):
                current = (-1)**m * comb(n, m) * (x**m / factorial(m))
                total += current
            return(total)

        def gamma(x, y):
            return ((2 * np.pi) * (((x - xc)**2 / a**2) + ((y - yc)**2 / b**2)))

        vfunc = np.vectorize(lgpoly)  # vectorize the lgpoly function

        g = gamma(x, y)
        lgp = vfunc(g)
        Clg = np.exp(-.5 * g) * lgp

        if norm is True:
            Clg = self._norm_channels(Clg)

        return Clg

    def _spatial_Gabor(self, b=2.5, theta=np.pi / 2, Lambda=20,
                       phi=np.pi / 4, gamma=1, k=np.pi, norm=True):
        '''
        Create 2D Gabor filter in spatial domain. 

        Code copied from: http://www.cs.rug.nl/~imaging/simplecell.html

        Parameters:
            b (type- float): specifies the spatial-frequency bandwidth of the
                filter. The bandwidth is specified in octaves.

            theta (type- float): orientation of the normal to the parallel
                stripes of a Gabor function. Specified in radians.

            Lambda (type- int): Wavelength of the sinusoidal (cosine) factor. 
                Specified in pixels. 

            phi (type- float): phase offset of cosine factor. Specified in 
                radians. 

            gamma (type- float):  spatial aspect ratio, specifies the 
                ellipticity of the support of the Gabor function. Specified as
                a ratio. 

        Return:

            Gabor filter of size "dimensions" defined in class __init__.

        '''
        sigma = Lambda * ((1 / k) * np.sqrt(np.log(2) / 2) *
                          (2**b + 1) / (2**b - 1))

        # import condition to check below...
        if (sigma > self.dim[0] / 6) or (sigma > self.dim[1] / 6):
            return None  # only return gabors that fit in the image
            # i.e. 2 standard deviations in each direction

        # Bounding Box
        x, y = self._euclidean2D()

        # Rotation
        x_theta = x * np.cos(theta) + y * np.sin(theta)
        y_theta = -x * np.sin(theta) + y * np.cos(theta)

        gaus = np.exp((-1 / (2 * sigma**2)) *
                      (x_theta**2 + y_theta**2 * gamma**2))
        sin = np.cos((2 * k) * (x_theta / Lambda) + phi)

        gabor = gaus * sin

        if norm is True:
            gabor = self._norm_channels(gabor)

        return gabor

    def get_channels(self):
        '''
        Create a filter bank for a specific filter type. 

        Parameters
        ----------
        ftype : TYPE, optional
            DESCRIPTION. The default is 'Gabor'.
        **kwargs : TYPE, dictionary
            DESCRIPTION. the dictionary maps all parameter value names for the 
            corresponding filter to sequences of those parameter values. This 
            creates a filter bank for the specific filter identified in ftype
            arguement above. 

        Returns
        -------
        Dictionary. mapping index of channel, in channel matrix defined in class
        CHO class below, to the channel filter. 

        '''
        # check whether a valid channel type has been chosen.
        checkChannelName(self.ftype)

        if len(self.channel_params.keys()) == 0:
            # default filter banks if no user defined parameters for filter are
            # specified in kwargs
            self.kwargs = defaultParams(
                self.ftype, self.channel_params["Gabor_PPD"])
        else:
            self.channel_params['LGauss_A'] = self.channel_params['LGauss_A'] if type(
                self.channel_params['LGauss_A']) is list else [self.channel_params['LGauss_A']]
            self.channel_params['LGauss_B'] = self.channel_params['LGauss_B'] if type(
                self.channel_params['LGauss_B']) is list else [self.channel_params['LGauss_B']]

        chnls = []  # initialize dictionary containing channels for CHO.
        # (start at 0) in channel matrix
        if self.ftype == 'DoG':
            for n in range(1, self.channel_params['DoG_N'] + 1):
                chnls.append(self._frequency_DoG(
                    n, alpha=self.channel_params['DoG_A'], q=self.channel_params['DoG_Q']))
        elif self.ftype == 'LGauss':
            for a in self.channel_params['LGauss_A']:
                for b in self.channel_params['LGauss_B']:
                    for n in range(1, self.channel_params['LGauss_N'] + 1):
                        chnls.append(
                            self._spatial_Laguerre_Gaussian(a=a, b=b, n=n))
        elif self.ftype == 'Gabor':
            ecc = 0
            self.channel_params["b"] = [1]
            self.channel_params["Phi"] = [0]
            self.channel_params["Gamma"] = [1]

            self.channel_params["Lambda"] = [(ecc + 1) * self.channel_params["Gabor_PPD"] / (
                2**i) for i in range(-1, self.channel_params["Gabor_O"])]  # convert to cycles per degree
            self.channel_params["Theta"] = deque(
                [(i / self.channel_params["Gabor_F"]) * np.pi for i in range(0, self.channel_params["Gabor_F"])])
            for b in self.channel_params['b']:
                for p in self.channel_params["Phi"]:
                    for g in self.channel_params["Gamma"]:
                        for la in self.channel_params["Lambda"]:
                            for t in self.channel_params['Theta']:
                                chnls.append(self._spatial_Gabor(
                                             b=b, theta=t, Lambda=la, phi=p,
                                             gamma=g, norm=True))
        else:
            pass
        return chnls
