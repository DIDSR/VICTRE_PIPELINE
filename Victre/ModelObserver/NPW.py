from . import ModelObserver
import numpy as np
from .LinSpace import LinSpace
import matplotlib.pyplot as plt
from .Exceptions import checkChannelName, checkChannelArgs, defaultParams
from collections import deque
from scipy.special import comb, factorial
import scipy


class NPW(ModelObserver):
    def __init__(self, signal=None, signal_present_samples=None, signal_absent_samples=None,
                 convolutional=False):
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
            self.dimensions = len(signal.shape)
            self.samples = None
            self.signal = signal
        else:
            raise Exception("You need to define the signal and background power spectrum \
                 or provide samples. ")

        self.build_template()
        self.template_statistics = self._calculate_template_statistics()

    def build_template(self):
        '''
        This method builds the Channelized Hotelling Template with the given
        channels and signal.

        Returns
        -------
        None.

        '''
        self.template = self.signal
        self.fft_template = np.fft.fftn(np.fft.fftshift(self.template))
