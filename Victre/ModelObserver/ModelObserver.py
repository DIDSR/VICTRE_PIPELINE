import numpy as np
from scipy.stats import norm
import glob
import progressbar
from termcolor import cprint
import os
import h5py
import copy
from sklearn.utils import shuffle
from sklearn.metrics import roc_curve, auc
from .. import Exceptions


class ModelObserver:
    """
    Initializes the ModelObserver generic class. This should be done on the children classes, not here.

    :param signal_present_samples: List of images (2D or 3D) for signal-present trials.
    :param signal_absent_samples: List of images (2D or 3D) for signal-absent trials.
    :param signal: Provide your own signal profile, otherwise, samples will be used to calculate the signal profile.
    :param results_folder: Path to the VICTRE results folder if trials are to be extracted from there.
    :param modality: Modality to use when path to VICTRE results folder is given.
    :param training_ratio: ratio of samples used to retrain.

    :returns: None
    """

    def __init__(self, signal_present_samples=None, signal_absent_samples=None, signal=None, results_folder=None, modality="dm", training_ratio=1, subtract_mean_background=False):
        self.template = None
        self.fft_template = None
        self.template_statistics = dict()
        self.samples = dict()
        self.dimensions = 2
        self.training_ratio = training_ratio
        self.responses = None
        self.responses_MR = None
        self.signal_present_samples = signal_present_samples
        self.signal_absent_samples = signal_absent_samples

        if results_folder is not None:
            # samples not provided, try to get them from the results folder
            self.read_from_VICTRE(
                results_folder, modality, subtract_mean_background=subtract_mean_background)
        else:
            # extract the signal from samples or use the provided signal
            self.extract_from_samples(
                signal_present_samples, signal_absent_samples)

        if signal is not None:
            self.dimensions = len(signal.shape)
            self.signal = signal

    def estimate_NPS(self, samples):
        """
        This method estimates the noise power spectrum from samples in any dimension.
        :param samples: List of samples.

        :returns: Noise power spectrum matrix.
        """
        mean = np.mean(samples, axis=0)
        # mean = np.mean(samples, 0)
        nps = np.zeros(self.signal.shape)
        for sample in samples:
            f = np.fft.fftn(sample - mean)
            nps += np.abs(f)**2

        return nps / len(samples)

    def run_readers(self, num_readers, LKE=True, function="sum"):
        """
        Run a specified number of readers, retraining them every iteration.

        :param num_readers: Number of readers to run.
        :param LKE: Indicate a Location-Known Exactly (LKE) task or a search task.
        :param function: For search tasks (`LKE = false`), specify the function used to
                combine the results of the convolution. Valid options are `sum` or `max`.

        :returns: 2-tuple with the statistics of each reader and their corresponding responses
                to each processed trial.
        """
        self.responses_MR = [None] * num_readers
        self.statistics_MR = [None] * num_readers
        for reader in range(num_readers):
            self.statistics_MR[reader], self.responses_MR[reader] = self.get_performance(
                function=function,
                LKE=LKE)
            self.retrain()

        return self.statistics_MR, self.responses_MR

    def save_iMRMC(self, filename, scores=None):
        """
        Save the results of the model observer as an iMRMC format file.

        :param filename: Output file name to save the result.
        :param scores: Alternatively, specify the scores of various modalities to be saved in
                      the file instead of the current model observer's results.

        :returns: None.
        """
        if self.responses_MR is None and scores is None:
            raise Exceptions.VictreError("Model needs to be run first!")

        if scores is None:
            scores = {}
            scores["modality"] = [sub['scores'] for sub in self.responses_MR]
            num_readers = len(self.responses_MR)
        else:
            self.responses_MR = []
            num_readers = len(scores[list(scores.keys())[0]])
        lines = []

        lines.append("Simulated iMRMC input")
        lines.append("")
        lines.append("NR: {:d}".format(num_readers))
        lines.append("N0: {:d}".format(len(self.signal_absent_samples)))
        lines.append("N1: {:d}".format(len(self.signal_present_samples)))
        lines.append("NM: {:d}".format(len(scores)))
        lines.append("")
        lines.append("BEGIN DATA:")

        for t in range(len(self.signal_present_samples)):
            lines.append("-1,P_{:d},0,1".format(t))
        for t in range(len(self.signal_absent_samples)):
            lines.append("-1,A_{:d},0,0".format(t))

        for modality in scores:
            for idx, responses in enumerate(scores[modality]):
                for t in range(len(responses["present"][0])):
                    lines.append("{:d},P_{:d},{:s},{:}".format(idx,
                                                               responses["present"][0][t],
                                                               modality,
                                                               responses["present"][1][t]))
                for t in range(len(responses["absent"][0])):
                    lines.append("{:d},A_{:d},{:s},{:}".format(idx,
                                                               responses["absent"][0][t],
                                                               modality,
                                                               responses["absent"][1][t]))

        with open(filename, "w") as outfile:
            outfile.write("\n".join(lines))

    def run(self, images=None, LKE=True, function="sum"):
        """
        Computes the result for the model observer template in a LKE or Search task on a set of images.

        :param images: List of images to be processed by the model observer.
        :param LKE: Indicate a Location-Known Exactly (LKE) task or a search task.
        :param function: For search tasks (`LKE = false`), specify the function used to
                combine the results of the convolution. Valid options are `sum` or `max`.

        :returns: List with the resulting statistic of the model observer computation for each image.
        """

        auto = False
        if images is None:
            auto = True
            images = self.testing_absent + self.testing_present

        elif type(images) is not list:
            images = [images]

        responses = np.zeros(len(images))
        LR = np.zeros(len(images))
        max_location = np.zeros((len(images),) + (len(images[0].shape),))
        max_location_LR = np.zeros((len(images),) + (len(images[0].shape),))

        for i in range(len(images)):
            image = images[i]
            if LKE:
                responses[i] = np.sum(image * self.template)
                LR[i] = responses[i]
            else:
                orig_shape = self.template.shape
                if len(self.template.shape) == 2:
                    image = np.pad(copy.deepcopy(image),
                                   ((self.template.shape[0], self.template.shape[0]),
                                    (self.template.shape[1], self.template.shape[1])),
                                   mode='constant')
                    if self.template.shape[0] != image.shape[0] or self.template.shape[1] != image.shape[1]:
                        self.template = np.fft.fftshift(np.pad(self.template, (
                            (np.ceil((image.shape[0] - self.template.shape[0]) / 2).astype(int),
                             np.floor((image.shape[0] - self.template.shape[0]) / 2).astype(int)),
                            (np.ceil((image.shape[1] - self.template.shape[1]) / 2).astype(int),
                             np.floor((image.shape[1] - self.template.shape[1]) / 2).astype(int))
                        ), mode="constant", constant_values=np.median(self.template)))
                        self.fft_template = np.fft.fftn(self.template)
                else:
                    image = np.pad(copy.deepcopy(image),
                                   ((self.template.shape[0], self.template.shape[0]),
                                    (self.template.shape[1],
                                     self.template.shape[1]),
                                    (self.template.shape[2], self.template.shape[2])),
                                   mode='constant')
                    if self.template.shape[0] != image.shape[0] or self.template.shape[1] != image.shape[1] or self.template.shape[2] != image.shape[2]:
                        self.template = np.fft.fftshift(np.pad(self.template, (
                            (np.ceil((image.shape[0] - self.template.shape[0]) / 2).astype(int),
                             np.floor((image.shape[0] - self.template.shape[0]) / 2).astype(int)),
                            (np.ceil((image.shape[1] - self.template.shape[1]) / 2).astype(int),
                             np.floor((image.shape[1] - self.template.shape[1]) / 2).astype(int)),
                            (np.ceil((image.shape[2] - self.template.shape[2]) / 2).astype(int),
                             np.floor((image.shape[2] - self.template.shape[2]) / 2).astype(int))
                        ), mode="constant", constant_values=np.median(self.template)))
                        self.fft_template = np.fft.fftn(self.template)

                fft_image = np.fft.fftn(image)
                template_responses = np.real(
                    np.fft.ifftn(fft_image * self.fft_template))

                if len(self.template.shape) == 2:
                    self.template_responses = template_responses[orig_shape[0]:template_responses.shape[0] - orig_shape[0],
                                                                 orig_shape[1]:template_responses.shape[1] - orig_shape[1]]
                else:
                    self.template_responses = template_responses[orig_shape[0]:template_responses.shape[0] - orig_shape[0],
                                                                 orig_shape[1]:template_responses.shape[1] - orig_shape[1],
                                                                 orig_shape[2]:template_responses.shape[2] - orig_shape[2]]

                if self.template_statistics is not None:
                    lhP = norm.pdf(self.template_responses,
                                   loc=self.template_statistics['mean_p'],
                                   scale=self.template_statistics['std'])
                    lhA = norm.pdf(self.template_responses,
                                   loc=self.template_statistics['mean_a'],
                                   scale=self.template_statistics['std'])

                    try:
                        with np.errstate(divide='ignore'):
                            self.LR = lhP / lhA
                            self.LR[self.LR == np.inf] = np.nanmax(
                                self.LR[np.isfinite(self.LR)])
                            self.LR[np.isnan(self.LR)] = np.nanmin(
                                self.LR)
                    except:
                        self.LR = lhP
                else:
                    self.LR = np.array([0])

                if (function == "sum"):
                    responses[i] = np.nansum(self.template_responses)
                    LR[i] = np.nansum(self.LR)
                else:
                    responses[i] = np.nanmax(self.template_responses)
                    max_location[i] = np.unravel_index(
                        np.argmax(self.template_responses), self.template_responses.shape)

                    LR[i] = np.nanmax(self.LR)
                    max_location_LR[i] = np.unravel_index(
                        np.argmax(self.LR), self.LR.shape)

        if auto:
            return responses[len(self.testing_absent):], responses[:len(self.testing_absent)]
        else:
            return responses, max_location, LR, max_location_LR

    def get_performance(self, LKE=True, function="sum"):
        """
        Calculates model observer's performance for the loaded samples.

        :param LKE: Indicate a Location-Known Exactly (LKE) task or a search task.
        :param function: For search tasks (`LKE = false`), specify the function used to
        combine the results of the convolution. Valid options are `sum` or `max`.
        """
        responses = {"present": [[], []], "absent": [[], []]}

        responses["present"][0] = self.idx_shuffled_testing_present
        responses["present"][1] = self.run(
            self.testing_present, function, LKE)[2]

        responses["absent"][0] = self.idx_shuffled_testing_absent
        responses["absent"][1] = self.run(
            self.testing_absent, function, LKE)[2]

        scores = np.append(np.array(responses["present"][1]),
                           np.array(responses["absent"][1]))
        truth = np.append([[1] * len(responses["present"][1])],
                          [[0] * len(responses["absent"][1])])
        fpr, tpr, _ = roc_curve(truth, scores)

        dP = (np.mean(responses["present"][1]) - np.mean(responses["absent"][1])) / (np.sqrt(
            0.5 * np.var(responses["present"][1]) + 0.5 * np.var(responses["absent"][1])))

        return {"fpr": fpr, "tpr": tpr, "auc": auc(fpr, tpr), "dP": dP}, responses

    def _calculate_template_statistics(self):
        """
        Calculates the analytical template statistics.

        :returns: Dictionary with the mean response for present samples,
                the mean response for absent samples and the standard deviation.
        """
        if self.samples is not None and self.samples["absent"].shape[-1] > 1 and self.samples["absent"].shape[-1] > 1:
            if self.dimensions == 2:
                sp = self.samples["present"].reshape(
                    (self.samples["present"].shape[0] * self.samples["present"].shape[1],
                     self.samples["present"].shape[2]))

                sa = self.samples["absent"].reshape(
                    (self.samples["absent"].shape[0] * self.samples["absent"].shape[1],
                     self.samples["absent"].shape[2]))

                template_1D = self.template.reshape(
                    (self.template.shape[0] * self.template.shape[1], 1))
            else:
                sp = self.samples["present"].reshape(
                    (self.samples["present"].shape[0] * self.samples["present"].shape[1] * self.samples["present"].shape[2],
                     self.samples["present"].shape[3]))

                sa = self.samples["absent"].reshape(
                    (self.samples["absent"].shape[0] * self.samples["absent"].shape[1] * self.samples["present"].shape[2],
                     self.samples["absent"].shape[3]))

                template_1D = self.template.reshape(
                    (self.template.shape[0] * self.template.shape[1] * self.template.shape[2], 1))

            t_present = np.dot(template_1D.transpose(), sp)
            t_absent = np.dot(template_1D.transpose(), sa)

            t_all = np.concatenate((t_present, t_absent), axis=1)

            ret = {"mean_p": np.mean(t_present),
                   "mean_a": np.mean(t_absent),
                   "std": np.std(t_all)}

        else:
            ret = None

        return ret

    @ staticmethod
    def _euclidean2D(dim):
        """
        Creates a linear space.

        :param dim: 2-tuple with the dimensions.

        :returns: 2-tuple with the linear space.
        """
        halfc = int(np.ceil((dim[1] - 1) / 2))
        halfr = int(np.ceil((dim[0] - 1) / 2))
        cols = np.linspace(-halfc, halfc - 1, dim[1])
        rows = np.linspace(halfr, -1 * (halfr - 1), dim[0])
        x, y = np.meshgrid(cols, rows)
        return (x, y)

    @ staticmethod
    def _euclidean3D(dim):
        """
        Creates a linear space in 3D.

        :param dim: 3-tuple with the dimensions.

        :returns: 3-tuple with the linear space.
        """
        halfc = int(np.ceil((dim[1] - 1) / 2))
        halfr = int(np.ceil((dim[0] - 1) / 2))
        halfd = int(np.ceil((dim[2] - 1) / 2))
        cols = np.linspace(-halfc, halfc - 1, dim[1])
        rows = np.linspace(halfr, -1 * (halfr - 1), dim[0])
        depth = np.linspace(halfd, -1 * (halfd - 1), dim[2])
        x, y, z = np.meshgrid(cols, rows, depth)
        return (x, y, z)

    def retrain(self, training_ratio=None):
        """
        Retrains the model from a different set of samples.

        :param training_ratio: ratio of samples used to retrain.
        :returns: None.
        """
        if training_ratio is not None:
            self.training_ratio = training_ratio

        self.extract_from_samples(
            self.signal_present_samples, self.signal_absent_samples)
        self.build_template()

    def extract_from_samples(self, signal_present_samples, signal_absent_samples):
        """
        Extracts the signal only from present and absent samples and saves it as class variables.

        :param signal_present_samples: List of signal-present images.
        :param signal_absent_samples: List of signal-absent images.

        :returns: None.
        """

        shuffled_present, idx_shuffled_present = shuffle(
            signal_present_samples, range(len(signal_present_samples)))
        shuffled_absent, idx_shuffled_absent = shuffle(
            signal_absent_samples, range(len(signal_absent_samples)))

        self.training_present = shuffled_present[:int(
            self.training_ratio * len(shuffled_present))]
        self.training_absent = shuffled_absent[:np.max((1, int(
            self.training_ratio * len(shuffled_absent))))]

        if self.training_ratio < 1 and self.training_ratio > 0:
            self.testing_present = shuffled_present[int(
                self.training_ratio * len(shuffled_present)):]
            self.testing_absent = shuffled_absent[int(
                self.training_ratio * len(shuffled_absent)):]

            self.idx_shuffled_testing_present = idx_shuffled_present[int(
                self.training_ratio * len(shuffled_present)):]
            self.idx_shuffled_testing_absent = idx_shuffled_absent[int(
                self.training_ratio * len(shuffled_absent)):]

        else:
            self.testing_present = self.training_present
            self.testing_absent = self.training_absent
            self.idx_shuffled_testing_present = idx_shuffled_present
            self.idx_shuffled_testing_absent = idx_shuffled_absent

        self.dimensions = len(signal_present_samples[0].shape)
        self.samples = {"absent": np.stack(self.training_absent, axis=self.dimensions),
                        "present": np.stack(self.training_present, axis=self.dimensions)}
        mean = {"absent": np.mean(self.samples["absent"], axis=self.dimensions),
                "present": np.mean(self.samples["present"], axis=self.dimensions)}
        self.signal = mean["present"] - mean["absent"]

        self.nps = self.estimate_NPS(signal_absent_samples)
        self.nps += self.estimate_NPS(signal_present_samples)
        self.nps /= 2

    def read_from_VICTRE(self, results_folder, modality, subtract_mean_background=True):
        """
        Reads the samples from a VICTRE results folder structure and extracts the signal.

        :param results_folder: Path to the results folder from the VICTRE pipeline.
        :param modality: Modality to load (dm or dbt)

        :returns: None.
        """
        cases = glob.glob("{:s}/*/".format(results_folder))

        self.signal_present_samples = []
        self.signal_absent_samples = []

        for n in progressbar.progressbar(range(len(cases))):
            try:
                seed = int(cases[n].split('/')[-2])
                if not os.path.isfile("{:s}/{:d}/ROIs.h5".format(results_folder, seed)):
                    cprint(
                        "\nROIs not found, skipping case {:d}".format(seed), 'red')
                    continue
                all_ROIs = h5py.File(
                    "{:s}/{:d}/ROIs.h5".format(results_folder, seed), "r")

                for idx, lesion_type in enumerate(all_ROIs["dm"]["lesion_type"]):
                    roi = all_ROIs[modality][str(idx)][()]

                    if subtract_mean_background:
                        roi = roi - np.abs(np.mean(roi))

                    # if modality == "dbt":  # DBT reconstruction is transposed
                    #     roi = np.transpose(roi, (0, 2, 1))

                    # only if roi is the same size we can append
                    if lesion_type < 0 and (len(self.signal_absent_samples) == 0 or len(self.signal_absent_samples) > 0 and self.signal_absent_samples[0].shape == roi.shape):
                        self.signal_absent_samples.append(roi)
                    elif lesion_type > 0 and (len(self.signal_present_samples) == 0 or len(self.signal_present_samples) > 0 and self.signal_present_samples[0].shape == roi.shape):
                        self.signal_present_samples.append(roi)
            except:
                cprint(
                    "\nError reading ROIs, skipping case {:d}".format(seed), 'red')

        if len(self.signal_absent_samples) == 0:
            self.signal_absent_samples = [
                np.zeros(self.signal_present_samples[0].shape)]

        self.extract_from_samples(
            self.signal_present_samples, self.signal_absent_samples)
