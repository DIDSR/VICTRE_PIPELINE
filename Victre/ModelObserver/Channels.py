from collections import deque
from scipy.special import comb, factorial
import numpy as np


class Channels:
    def __init__(self, space, channel_params, channel_type):
        """!
        Initialize the channels for the CHO model.

        @param space 2-tuple with the space in which the channels will be generated (e.g. `self._euclidean2D(SIZE)`).
        @param channel_params The dictionary maps all parameter value names for the 
            corresponding channel to sequences of those parameter values. This 
            creates a channel set for the specified channel_type arguement below. 
        @param channel_type Type of the channels to be created for the CHO: Gabor, LGauss or DoG.
        @return None

        """
        self.space = space
        self.channel_type = channel_type
        self.channel_params = channel_params

    def _norm_channels(self, channel):
        """!
        Normalize the values in the 2D channel template. 

        @param channel Channel to be normalized
        @return Normalized channel

        """
        return (channel / np.sum(channel))

    def _frequency_DoG(self, q=1.67, alpha=1.4, s0=.005, n=1, norm=True):
        """!
        Create a Difference of Gaussians (DoG) channel in the frequency domain

        @param q Parameter q for the DoG formula.
        @param alpha Parameter alpha for the DoG formula.
        @param s0 Parameter s0 for the DoG formula.
        @param n Channel number.
        @param norm Activate channel normalization.
        @return DoG channel in the frequency domain with the DC component in the center.

        """
        u, v = self.space[0] / self.dim[1], self.space[1] / self.dim[0]
        freqs = np.sqrt(u**2 + v**2)  # 2D spatial frequency domain

        sj = s0 * (alpha**n)
        g1 = np.exp((-1 / 2) * (freqs / (q * sj))**2)
        g2 = np.exp((-1 / 2) * (freqs / sj)**2)

        dog = g1 - g2
        if norm is True:
            dog = self._norm_channels(dog)
        return dog

    def _spatial_Laguerre_Gaussian(self, a=30, b=25, c=30, n=5, norm=True, signal=None):
        """!
        Generate a Laguerre—Gauss (LGauss) channel in the spatial domain. 

        @param q Parameter a for the LGauss formula.
        @param b Parameter b for the LGauss formula.
        @param n Channel number.
        @param norm Activate channel normalization.
        @return LGauss channel in the spatial domain.
        """
        xc, yc, zc = 0, 0, 0

        def lgpoly(x):
            total = 0
            for m in range(n):
                current = (-1)**m * comb(n, m) * (x**m / factorial(m))
                total += current
            return(total)

        # def gamma(x, y):
        #     return ((2 * np.pi) * (((x - xc)**2 / a**2) + ((y - yc)**2 / b**2)))

        def gamma(x, y, z=0):
            return ((2 * np.pi) * (((x - xc)**2 / a**2) + ((y - yc)**2 / b**2) + ((z - zc)**2 / c**2)))

        vfunc = np.vectorize(lgpoly)  # vectorize the lgpoly function

        if signal is not None:
            if np.max(signal) != np.min(signal):
                img = (1 * (signal - np.min(signal)) /
                       (np.max(signal) - np.min(signal)))
                img = 1 - img
            else:
                img = signal
            if len(self.space) == 3:
                linear = self.space[0] * \
                    img, self.space[1] * img, self.space[2] * img
            else:
                linear = self.space[0] * img, self.space[1] * img
        else:
            if len(self.space) == 3:
                linear = self.space[0], self.space[1], self.space[2]
            else:
                linear = self.space[0], self.space[1]

        if len(self.space) == 3:
            g = gamma(linear[0], linear[1], linear[2])
        else:
            g = gamma(linear[0], linear[1])

        lgp = vfunc(g)
        Clg = np.exp(-.5 * g) * lgp

        if norm is True:
            Clg = self._norm_channels(Clg)

        return Clg

    def _spatial_Gabor(self, b=2.5, theta=3.14159265359 / 2, lmbd=20,
                       phi=3.14159265359 / 4, gamma=1, k=3.14159265359, norm=True):
        """!
        Generate a Gabor channel in the spatial domain. 

        @param b Specifies the spatial-frequency bandwidth of the
                filter. The bandwidth is specified in octaves.
        @param theta Orientation of the normal to the parallel
                stripes of a Gabor function. Specified in radians.
        @param lmbd Lambda. Wavelength of the sinusoidal (cosine) factor. 
                Specified in pixels.
        @param phi Phase offset of cosine factor. Specified in 
                radians. 
        @param gamma Spatial aspect ratio, specifies the ellipticity 
                of the support of the Gabor function. Specified as
                a ratio.
        @param k Parameter k for the Gabor formula.
        @param norm Activate channel normalization.

        @return Gabor channel in the spatial domain.
        """

        sigma = lmbd * ((1 / k) * np.sqrt(np.log(2) / 2) *
                        (2**b + 1) / (2**b - 1))

        # import condition to check below...
        if (sigma > self.space[0].shape[0] / 6) or (sigma > self.space[0].shape[1] / 6):
            return None  # only return gabors that fit in the image
            # i.e. 2 standard deviations in each direction

        # Rotation
        x_theta = self.space[0] * np.cos(theta) + self.space[1] * np.sin(theta)
        y_theta = -self.space[0] * \
            np.sin(theta) + self.space[1] * np.cos(theta)

        gaus = np.exp((-1 / (2 * sigma**2)) *
                      (x_theta**2 + y_theta**2 * gamma**2))
        sin = np.cos((2 * k) * (x_theta / lmbd) + phi)

        gabor = gaus * sin

        if norm is True:
            gabor = self._norm_channels(gabor)

        return gabor

    @staticmethod
    def checkChannelName(name):
        names = ['DoG', 'Gabor', 'LGauss', 'LGauss_signal']
        if name not in names:
            er = """This filter type, -> "{0}" <- has not been implemented yet.
            Acceptable arguments for channel type are: '{1}', '{2}', '{3}' and '{4}'.
            """.format(name, *names)
            raise Exception(er)
        else:
            pass

    def get_channels(self):
        """!
        Creates the channel set. 

        @return List of channels.
        """
        # check whether a valid channel type has been chosen.
        Channels.checkChannelName(self.channel_type)

        if self.channel_params is None or len(self.channel_params.keys()) == 0:
            # default channels if no user defined parameters
            self.channel_params = Channels.defaultParams(self.channel_type)
        elif "LGauss" in self.channel_type:
            self.channel_params['LGauss_A'] = self.channel_params['LGauss_A'] if type(
                self.channel_params['LGauss_A']) is list else [self.channel_params['LGauss_A']]
            self.channel_params['LGauss_B'] = self.channel_params['LGauss_B'] if type(
                self.channel_params['LGauss_B']) is list else [self.channel_params['LGauss_B']]
            if "LGauss_C" in self.channel_params:
                self.channel_params['LGauss_C'] = self.channel_params['LGauss_C'] if type(
                    self.channel_params['LGauss_C']) is list else [self.channel_params['LGauss_C']]

        chnls = []  # initialize dictionary containing channels for CHO.

        if self.channel_type == 'DoG':
            for n in range(1, self.channel_params['DoG_N'] + 1):
                chnls.append(self._frequency_DoG(
                    alpha=self.channel_params['DoG_A'],
                    q=self.channel_params['DoG_Q'],
                    n=n))
        elif "LGauss" in self.channel_type:
            if "LGauss_signal" not in self.channel_params:
                self.channel_params['LGauss_signal'] = None
            if "LGauss_C" not in self.channel_params:
                self.channel_params["LGauss_C"] = [1]
            for a in self.channel_params['LGauss_A']:
                for b in self.channel_params['LGauss_B']:
                    for c in self.channel_params["LGauss_C"]:
                        for n in range(1, self.channel_params['LGauss_N'] + 1):
                            chnls.append(
                                self._spatial_Laguerre_Gaussian(a=a, b=b, c=c, n=n, signal=self.channel_params['LGauss_signal']))
        elif self.channel_type == 'Gabor':
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

        return chnls

    @staticmethod
    def defaultParams(channel_type):
        # default to 45 pixels per degree visual angle for Gabor channels, see "Lambda" param below.
        pixels_per_degree = 45
        ecc = 1  # this would be the eccentricity for Gabor channels
        # convert wavelength of sinusoid factor to cycles per pixel and scale by c
        default_parameters = {
            'DoG': {
                'N': 10,
            },
            'Gabor': {
                'b': [1],
                'Theta': deque([(i / 8) * np.pi for i in range(0, 8)]),
                # convert to cycles per degree
                'Lambda': [ecc * pixels_per_degree / (2**i) for i in range(-1, 5)],
                'Phi': [0],
                'Gamma': [1],
            },
            'LGauss': {
                'LGauss_A': [30],
                'LGauss_B': [30],
                'LGauss_N': 5,
            },
            'LGauss_signal': None
        }
        return(default_parameters[channel_type])
