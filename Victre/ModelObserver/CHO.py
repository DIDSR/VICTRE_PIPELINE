from . import ModelObserver
import numpy as np
from .Channels import Channels
import scipy
import copy


class CHO(ModelObserver):

    def __init__(self, channel_type=None, channels=None, signal=None, nps=None, signal_present_samples=None,
                 signal_absent_samples=None, results_folder=None, modality="dm", training_ratio=1, channel_filter=None, channel_params=None, subtract_mean_background=False):
        """!
        This method initializes a Channelized Hotelling Observer (CHO) given a specific
        channel type. It can receive the signal or a set of training present and
        absent images.

        @param channel_type Type of the channels to be created for the CHO: Gabor, LGauss or DoG.
        @param channels Alternatively, you can provide your own channels.
        @param signal Provide your own signal profile, otherwise, samples will be used to calculate the signal profile.
        @param nps Noise power spectrum used to filter the channels.
        @param signal_present_samples List of images (2D or 3D) for signal-present trials.
        @param signal_absent_samples List of images (2D or 3D) for signal-absent trials.
        @param results_folder Path to the VICTRE results folder if trials are to be extracted from there.
        @param modality Modality to use when path to VICTRE results folder is given.
        @param training_ratio ratio of samples used to retrain.
        @param channel_filter Use of filtered channels.
        @param channel_params Dictionary with the parameters for the channel generation.

        @return None
        """

        super().__init__(signal_present_samples,
                         signal_absent_samples, signal, results_folder, modality, training_ratio, subtract_mean_background)

        self.channel_params = channel_params

        self.signal = self.signal / \
            np.abs(np.sum(self.signal))  # normalize signal

        if nps is not None:
            self.samples = None
            self.nps = nps

        self.fft_signal = np.fft.fftn(self.signal)

        if channel_filter is not None and "3D" in channel_filter and len(self.signal.shape) == 3:
            space = self._euclidean3D(self.signal.shape)
        else:
            space = self._euclidean2D(self.signal.shape[-2:])

        if channel_filter is not None and "LGauss" in channel_filter:
            self.channel_params["LGauss_signal"] = self.signal

        # get channels ready for the CHO
        if channels is None and channel_type is not None:
            # c = (1 + a * (ecci)**b)
            self.channels = Channels(space,
                                     self.channel_params,
                                     channel_type).get_channels()
            self.channels = [ch for ch in self.channels if ch is not None]
        else:
            self.channels = channels

        # repeat channels in 3D if they are not already 3D
        if self.dimensions == 3 and len(self.channels[0].shape) != 3 and channel_filter != "convolutional":
            if channel_filter is not None and "LGauss" in channel_filter:
                for ch in range(len(self.channels)):
                    self.channels[ch] = np.zeros(self.signal.shape)

                for slice in range(self.signal.shape[0]):
                    self.channel_params["LGauss_signal"] = self.signal[slice, :, :]

                    channels = Channels(space,
                                        self.channel_params,
                                        channel_type).get_channels()
                    channels = [ch for ch in channels if ch is not None]

                    for ch in range(len(self.channels)):
                        self.channels[ch][slice, :, :] = channels[ch]
            else:
                for ch in range(len(self.channels)):
                    # transpose channels if modality is DBT
                    self.channels[ch] = np.dstack(
                        [self.channels[ch]] * self.signal.shape[0]).transpose()

        if channel_filter is not None and channel_filter == "convolutional":
            for ch in range(len(self.channels)):
                this_channel = copy.deepcopy(self.channels[ch])**2
                self.channels[ch] = np.zeros(self.signal.shape)
                if self.dimensions == 3:
                    for slice in range(self.signal.shape[0]):
                        self.channels[ch][slice, :, :] = scipy.ndimage.convolve(this_channel,
                                                                                self.signal[slice, :, :],
                                                                                mode='nearest')
                else:
                    self.channels[ch] = scipy.ndimage.convolve(this_channel, self.signal,
                                                               mode='nearest')
                # else:
                #     self.channels[ch] = np.real(np.fft.ifftn(np.multiply(
                #         np.abs(
                #             np.fft.fftn(np.fft.fftshift(self.channels[ch]))) ** 2,
                #         self.fft_signal
                #     )))

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
        """!
        This method builds the Channelized Hotelling template with the given
        channels and signal.

        @return None

        """
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
            if np.all(v_absent == 0):
                self.K = np.cov(v_present)
            else:
                self.K = (np.cov(v_present) + np.cov(v_absent)) / 2

        else:
            # if we don't have samples, we should have the signal by itself
            # and the covariance of the background
            self.s = np.dot(self.channels_1D.transpose(), self.signal)
            self.K = self._analytic_covariance_matrix()

        # self.K += np.eye(self.K.shape[0]) * self.K

        # linear weights for the Hotelling template
        self.weights = np.dot(self.s, np.linalg.pinv(self.K))

        # synthesize a spatial template or kernel
        self.template = np.reshape(np.dot(self.weights, self.channels_1D.transpose()),
                                   self.signal.shape)

        # eliminate edge artifacts
        for n in range(10):
            if len(self.template.shape) == 3:
                # self.template[n,:,:]=self.template[n,:,:]*n/10
                # self.template[-n,:,:]=self.template[-n,:,:]*n/10
                self.template[:, n, :] = self.template[:, n, :] * n / 10
                self.template[:, -n, :] = self.template[:, -n, :] * n / 10
                self.template[:, :, n] = self.template[:, :, n] * n / 10
                self.template[:, :, -n] = self.template[:, :, -n] * n / 10
            else:
                self.template[n, :] = self.template[n, :] * n / 10
                self.template[-n, :] = self.template[-n, :] * n / 10
                self.template[:, n] = self.template[:, n] * n / 10
                self.template[:, -n] = self.template[:, -n] * n / 10

        # self.template = (self.template - self.template.min())
        with np.errstate(divide='raise', invalid='raise'):
            self.template /= np.abs(self.template.sum())
        # we also save it in the frequency domain
        self.fft_template = np.fft.fftn(np.fft.fftshift(self.template))

    def _analytic_covariance_matrix(self):
        """!
        This method calculates the channel covariance matrix given
        a noise power spectrum for the background noise.

        @return 2D array containing the covariance matrix for the channels.

        """
        K = np.zeros((len(self.channels), len(self.channels)))
        for i in range(len(self.channels)):
            chFiltered = np.fft.fft2(self.channels[i]) * self.nps
            chSpatialD = np.fft.ifft2(chFiltered).real.flatten().reshape(1, -1)
            for j in range(len(self.channels)):
                ch2 = self.channels[j].flatten().reshape(1, -1)
                K[i, j] = chSpatialD @ ch2.transpose()
                K[j, i] = K[i, j]
        return K
