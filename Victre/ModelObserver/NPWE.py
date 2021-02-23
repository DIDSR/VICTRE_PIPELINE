from . import ModelObserver
import numpy as np
import matplotlib.pyplot as plt
from .Exceptions import checkChannelName, checkChannelArgs, defaultParams
from collections import deque
from scipy.special import comb, factorial
import scipy


class NPWE(ModelObserver):
    def __init__(self, signal=None, signal_present_samples=None, signal_absent_samples=None,
                 convolutional=False, pixels_per_degree=50, alpha=1.9, beta=2, gamma=0.5):
        '''
        This method initializes a Channelized Hotelling Observer (CHO) given a specific
        channel type. It can receive the signal or a set of training present and
        absent images.

        Parameters
        ----------
        signal : 2D or 3D array
            You can provide your own signal profile.
        signal_present_samples : list
            List of images (2D or 3D) for signal-present trials.
        signal_absent_samples : list
            List of images (2D or 3D) for signal-absent trials.
        convolutional : string or None
            Set to None to use FFT convolution, set to "convolve" to use spatial-kernel
            convolution.
        pixels_per_degree : float
            Pixels per degree of visual angle to generate the eye filter.

        Returns
        -------
        None.

        '''
        super().__init__()

        # extract the signal from sasmples or use the provided signal
        if signal_present_samples is not None:
            self.dimensions = len(signal_present_samples[0].shape)
            self.samples = {"absent": np.stack(signal_absent_samples, axis=self.dimensions),
                            "present": np.stack(signal_present_samples, axis=self.dimensions)}
            mean = {"absent": np.mean(self.samples["absent"], axis=self.dimensions),
                    "present": np.mean(self.samples["present"], axis=self.dimensions)}
            self.signal = mean["present"] - mean["absent"]
        elif signal is not None:

            self.samples = None
            self.signal = signal
        else:
            raise Exception("You need to define the signal and background power spectrum \
                 or provide samples. ")

        self.pixels_per_degree = pixels_per_degree
        self.eye_filter_parameters = {
            "alpha": alpha, "beta": beta, "gamma": gamma}

        self.build_template()
        self.template_statistics = self._calculate_template_statistics()

    def eye_filter(self, r):
        return (r ** self.eye_filter_parameters["alpha"] *
                np.exp(-self.eye_filter_parameters["beta"] * r ** self.eye_filter_parameters["gamma"])) ** 2

    def build_template(self):
        '''
        This method builds the Channelized Hotelling Template with the given
        channels and signal.

        Returns
        -------
        None.

        '''
        x, y = self._euclidean2D()
        u, v = x / self.dim[1], y / self.dim[0]
        freqspace = np.sqrt(u**2 + v**2)  # 2D spatial frequency domain

        self.freqs = np.fft.fftshift(freqspace) * self.pixels_per_degree

        self.E = np.array(list(map(self.eye_filter, self.freqs)))

        if self.dimensions == 3:  # repeat channels in 3D
            self.E = np.dstack([self.E] * self.signal.shape[-1])

        self.template = np.real(np.fft.ifftn(
            np.fft.fftn(self.signal) * self.E))
        self.fft_template = np.fft.fftn(np.fft.fftshift(self.template))
