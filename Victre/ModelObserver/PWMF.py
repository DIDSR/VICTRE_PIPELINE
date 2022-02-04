from . import ModelObserver
import numpy as np


class PWMF(ModelObserver):
    def __init__(self, signal=None, nps=None, signal_present_samples=None, signal_absent_samples=None, results_folder=None, modality="dm", training_ratio=1, subtract_mean_background=False):
        """!
        This method initializes a PreWhitening Matched Filter (PWMF) model observer. It can receive 
        the signal or a set of training present and absent images.

        @param signal Provide your own signal profile, otherwise, samples will be used to calculate the signal profile.
        @param nps Noise power spectrum used to filter the channels. If not provided it will be estimated from trials.
        @param signal_present_samples List of images (2D or 3D) for signal-present trials.
        @param signal_absent_samples List of images (2D or 3D) for signal-absent trials.
        @param results_folder Path to the VICTRE results folder if trials are to be extracted from there.
        @param modality Modality to use when path to VICTRE results folder is given.
        @param training_ratio ratio of samples used to retrain.

        @return None
        """
        super().__init__(signal_present_samples,
                         signal_absent_samples, signal, results_folder, modality, training_ratio, subtract_mean_background)

        if nps is not None:
            self.samples = None
            self.nps = nps

        self.build_template()
        self.template_statistics = self._calculate_template_statistics()

    def build_template(self):
        """!
        This method builds the PWMF template.

        @return None

        """
        self.fft_template = np.fft.fftn(
            np.fft.fftshift(self.signal)) / self.nps  # / np.sqrt(self.nps)

        self.template = np.fft.ifftshift(
            np.real(np.fft.ifftn(self.fft_template)))
        # self.template = (self.template - self.template.min()) / (self.template.max()-self.template.min())
        self.template = self.template / np.abs(self.template.sum())
        self.fft_template = np.fft.fftn(np.fft.fftshift(self.template))
