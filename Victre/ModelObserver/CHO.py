from . import ModelObserver
from .LinSpace import LinSpace
import numpy as np
import matplotlib.pyplot as plt
from .Exceptions import checkChannelName, checkChannelArgs, defaultParams
from .Channels import Channels
from collections import deque
from scipy.special import comb, factorial
import scipy


class CHO(ModelObserver):

    def __init__(self, channel_type=None, channels=None, signal=None, nps=None, signal_present_samples=None,
                 signal_absent_samples=None, convolutional=False, channel_params=None):
        '''
        This method initializes a Channelized Hotelling Observer (CHO) given a specific
        channel type. It can receive the signal or a set of training present and
        absent images.

        Parameters
        ----------
        channel_type : string
            Type of the channels to be created for the CHO: Gabor, LGauss or DoG.
        channels : list
            Alternatively, you can provide your own channels.
        signal : 2D or 3D array
            You can provide your own signal profile.
        nps : 2D or 3D array
            Noise power spectrum used to generate the background noise.
        signal_present_samples : list
            List of images (2D or 3D) for signal-present trials.
        signal_absent_samples : list
            List of images (2D or 3D) for signal-absent trials.
        convolutional : string or None
            Set to None to use FFT convolution, set to "convolve" to use spatial-kernel
            convolution.
        pixels_per_degree : float
            Pixels per degree of visual angle to generate the channels.

        Returns
        -------
        None.

        '''
        super().__init__()

        self.channel_params = channel_params

        # extract the signal from samples or use the provided signal
        if signal_present_samples is not None:
            self.dimensions = len(signal_present_samples[0].shape)
            self.samples = {"absent": np.stack(signal_absent_samples, axis=self.dimensions),
                            "present": np.stack(signal_present_samples, axis=self.dimensions)}
            mean = {"absent": np.mean(self.samples["absent"], axis=self.dimensions),
                    "present": np.mean(self.samples["present"], axis=self.dimensions)}
            self.signal = mean["present"] - mean["absent"]

        if signal is not None:
            self.dimensions = len(signal.shape)
            self.signal = signal

        self.signal = self.signal / np.sum(self.signal)  # normalize signal

        if nps is not None:
            self.samples = None
            self.nps = nps

        self.fft_signal = np.fft.fftn(self.signal)

        # get channels ready for the CHO
        if channels is None and channel_type is not None:
            # c = (1 + a * (ecci)**b)
            self.channels = Channels(signal_present_samples[0].shape[-2:],
                                     self.channel_params, channel_type).get_channels()
            self.channels = [ch for ch in self.channels if ch is not None]
        else:
            self.channels = channels

        if self.dimensions == 3:  # repeat channels in 3D
            for ch in range(len(self.channels)):
                # transpose channels if modality is DBT
                self.channels[ch] = np.dstack(
                    [self.channels[ch]] * self.signal.shape[0]).transpose()

        if convolutional:
            for ch in range(len(self.channels)):
                if convolutional == "convolve":
                    self.channels[ch] = scipy.ndimage.filters.convolve(self.channels[ch],
                                                                       self.signal,
                                                                       mode='nearest')
                else:
                    self.channels[ch] = np.real(np.fft.ifftn(np.multiply(
                        np.abs(
                            np.fft.fftn(np.fft.fftshift(self.channels[ch]))) ** 2,
                        self.fft_signal
                    )))

        self.channels_1D = np.stack(self.channels, axis=len(self.signal.shape))

        if self.dimensions == 2:
            self.channels_1D = self.channels_1D.reshape(
                (self.channels_1D.shape[0] * self.channels_1D.shape[1],
                 self.channels_1D.shape[2]))
        else:
            self.channels_1D = self.channels_1D.reshape(
                (self.channels_1D.shape[0] * self.channels_1D.shape[1] * self.channels_1D.shape[2],
                 self.channels_1D.shape[3]))

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
        # if we have samples, we can do it this way
        if self.samples is not None:
            if self.dimensions == 2:
                sp = self.samples["present"].reshape(
                    (self.samples["present"].shape[0] * self.samples["present"].shape[1],
                     self.samples["present"].shape[2]))

                sa = self.samples["absent"].reshape(
                    (self.samples["absent"].shape[0] * self.samples["absent"].shape[1],
                     self.samples["absent"].shape[2]))
            else:
                sp = self.samples["present"].reshape(
                    (self.samples["present"].shape[0] * self.samples["present"].shape[1] * self.samples["present"].shape[2],
                     self.samples["present"].shape[3]))

                sa = self.samples["absent"].reshape(
                    (self.samples["absent"].shape[0] * self.samples["absent"].shape[1] * self.samples["present"].shape[2],
                     self.samples["absent"].shape[3]))

            # response of the channels to signal-present trials
            v_present = np.dot(self.channels_1D.transpose(), sp)

            # response of the channels to signal-absent trials
            v_absent = np.dot(self.channels_1D.transpose(), sa)

            # response of the channels to the signal
            self.s = np.mean(v_present, axis=1) - np.mean(v_absent, axis=1)

            # covariance matrix
            self.K = (np.cov(v_present) + np.cov(v_absent)) / 2
        else:
            # if we don't have samples, we should have the signal by itself
            # and the covariance of the background
            self.s = np.dot(self.channels_1D.transpose(), self.signal)
            self.K = self._analytic_covariance_matrix()

        # linear weights for the Hotelling template
        self.weights = np.dot(self.s, np.linalg.pinv(self.K))

        # synthesize a spatial template or kernel
        self.template = np.reshape(np.dot(self.weights, self.channels_1D.transpose()),
                                   self.signal.shape)

        # we also save it in the frequency domain
        self.fft_template = np.fft.fftn(np.fft.fftshift(self.template))

    def _analytic_covariance_matrix(self):
        '''
        This method calculates the channel covariance matrix given
        a noise power spectrum for the background noise.

        Returns
        -------
        K : array
            2D array containing the covariance matrix for the channels.

        '''
        K = np.zeros(2 * [len(self.channels)])
        for i in range(len(self.channels)):
            chFiltered = np.fft.fft2(self.channels[i]) * self.nps
            chSpatialD = np.fft.ifft2(chFiltered).real.flatten().reshape(1, -1)
            for j in range(len(self.channels)):
                ch2 = self.channels[j].flatten().reshape(1, -1)
                K[i, j] = chSpatialD @ ch2
                K[j, i] = K[i, j]
        return K
