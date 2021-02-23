import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import plotly.express as px


class ModelObserver:
    def __init__(self):
        '''
        Initialize class.

        Parameters
        ----------
        dimensions : TYPE, tuple containing 2 (or 3) integers,
        i.e. len(dimensions) == 2 (or 3).
            DESCRIPTION. (# rows, # columns, optional: # of slices). This
                          defines the space in which the model observer class
                          will operate in.

        Returns
        -------
        None.

        '''
        self.template = None
        self.fft_template = None
        self.template_statistics = dict()
        self.samples = dict()
        self.dimensions = 2

    def covariance_matrix(self):
        pass

    def run(self, images, function="sum", LKE=True, location=None):
        '''
        This method takes a list of images and computes the result for the model 
        observer template in a LKE or Search task.

        Parameters
        ----------
        images : list
            list of images to be seen by the model observer.
        function : string, optional
            string with the function to use when LKE = False, max or sum.
        LKE : boolean, optional
            indicates a LKE task (True) or Search task (False).
        location : tuple, optional
            2 (or 3) integers with the signal location in the image for a LKE task.

        Returns
        -------
        responses : list
            List with the result of the model observer computation for each image.
        '''

        # images = np.stack(images, axis=len(images[0].shape)-1)

        if type(images) is not list:
            images = [images]

        responses = np.zeros(len(images))
        for i in range(len(images)):
            image = images[i]
            if LKE:
                responses[i] = np.sum(image * self.template)
            else:
                fft_image = np.fft.fftn(image)
                template_responses = np.real(
                    np.fft.ifftn(fft_image * self.fft_template))

                if self.template_statistics is not None:
                    lhP = norm.pdf(template_responses,
                                   loc=self.template_statistics['mean_p'],
                                   scale=self.template_statistics['std'])
                    lhA = norm.pdf(template_responses,
                                   loc=self.template_statistics['mean_a'],
                                   scale=self.template_statistics['std'])

                    with np.errstate(divide='ignore'):
                        self.template_responses = lhP / lhA
                if (function == "sum"):
                    responses[i] = np.sum(template_responses)
                else:
                    responses[i] = np.max(template_responses)

        return responses

    def _calculate_template_statistics(self):
        if self.samples is not None:
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

            ret = {"mean_p": np.mean(t_present), "mean_a": np.mean(
                t_absent), "std": np.std(t_all)}

        else:
            ret = None

        return ret

    def _euclidean2D(self):
        halfc = math.ceil((self.dim[1] - 1) / 2)
        halfr = math.ceil((self.dim[0] - 1) / 2)
        cols = np.linspace(-halfc, halfc - 1, self.dim[1])
        rows = np.linspace(halfr, -1 * (halfr - 1), self.dim[0])
        x, y = np.meshgrid(cols, rows)
        return (x, y)
