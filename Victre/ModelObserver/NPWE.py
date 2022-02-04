from . import ModelObserver
import numpy as np


class NPWE(ModelObserver):
    def __init__(self, signal=None, signal_present_samples=None, signal_absent_samples=None, results_folder=None, modality="dm", training_ratio=1,
                 pixels_per_degree=50, alpha=1.9, beta=2, gamma=0.5):
        """!
        This method initializes a Non PreWhitening with Eye filter (NPWE) model observer. It can receive 
        the signal or a set of training present and absent images.

        @param signal Provide your own signal profile, otherwise, samples will be used to calculate the signal profile.
        @param signal_present_samples List of images (2D or 3D) for signal-present trials.
        @param signal_absent_samples List of images (2D or 3D) for signal-absent trials.
        @param results_folder Path to the VICTRE results folder if trials are to be extracted from there.
        @param modality Modality to use when path to VICTRE results folder is given.
        @param training_ratio ratio of samples used to retrain.
        @param pixels_per_degree Pixels per degree of visual angle to generate the eye filter.
        @param alpha Parameter alpha for the eye filter.
        @param beta Parameter beta for the eye filter.
        @param gamma Parameter gamma for the eye filter.

        @return None
        """

        super().__init__(signal_present_samples,
                         signal_absent_samples, signal, results_folder, modality, training_ratio)

        self.pixels_per_degree = pixels_per_degree
        self.eye_filter_parameters = {
            "alpha": alpha, "beta": beta, "gamma": gamma}

        self.build_template()
        self.template_statistics = self._calculate_template_statistics()

    def eye_filter(self, r):
        """!
        Generates the corresponding eye filter at distance `r`.

        @param r Radial distance to generate the eye filter.

        @return Eye filter response at distance `r`.
        """
        return (r ** self.eye_filter_parameters["alpha"] *
                np.exp(-self.eye_filter_parameters["beta"] * r ** self.eye_filter_parameters["gamma"])) ** 2

    def build_template(self):
        """!
        This method builds the NPWE template.

        @return None

        """
        x, y = self._euclidean2D(self.signal.shape[-2:])
        u, v = x / self.dim[1], y / self.dim[0]
        freqspace = np.sqrt(u**2 + v**2)  # 2D spatial frequency domain

        self.freqs = np.fft.fftshift(freqspace) * self.pixels_per_degree

        self.E = np.array(list(map(self.eye_filter, self.freqs)))

        if self.dimensions == 3:  # repeat channels in 3D
            self.E = np.dstack([self.E] * self.signal.shape[-1])

        self.template = np.real(np.fft.ifftn(
            np.fft.fftn(self.signal) * self.E))
        self.template = (self.template - self.template.min()) / \
            (self.template.max() - self.template.min())
        self.template = self.template / np.abs(self.template.sum())
        self.fft_template = np.fft.fftn(np.fft.fftshift(self.template))
