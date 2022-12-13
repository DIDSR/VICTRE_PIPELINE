import numpy as np
from scipy import ndimage
from Victre.ModelObserver import CHO, NPW, PWMF
import unittest


class TestModelObserver(unittest.TestCase):

    size = [128, 128, 128]
    contrast = 0.05
    present, absent = [], []
    s = None

    truth = {
        "NPW": [0.9228571428571428, 0.6648979591836734, 0.6540816326530613, 0.606938775510204, 0.7306122448979593,
                0.7453061224489795, 0.6461224489795918, 0.6691836734693877, 0.676734693877551, 0.7012244897959183],
        "PWMF": [1.0, 0.769795918367347, 0.7822448979591837, 0.7773469387755102, 0.8340816326530612,
                 0.7773469387755102, 0.8238775510204082, 0.7406122448979592, 0.753469387755102, 0.7559183673469386],
        "CHO_LG": [0.9906122448979592, 0.9916326530612245, 0.9812244897959184, 0.9891836734693877, 0.9867346938775511,
                   0.9924489795918366, 0.9593877551020409, 0.9861224489795918, 0.9834693877551021, 0.986530612244898],
        "CHO_DoG": [0.9989795918367348, 0.9985714285714286, 0.9993877551020408, 1.0, 0.996734693877551,
                    0.9997959183673469, 0.9989795918367348, 0.9997959183673469, 0.9997959183673469, 0.9975510204081632],
        "CHO_Gabor": [0.9789795918367348, 0.9936734693877551, 0.9857142857142858, 0.9871428571428572, 0.9897959183673469,
                      0.9836734693877551, 0.993673469387755, 0.9918367346938777, 0.9914285714285714, 0.9987755102040816]
    }

    @classmethod
    def setUpClass(cls):
        np.random.seed(7127433)
        for _ in range(100):
            cls.absent.append(cls._genNoiseNPS(cls.size[:2], 2, 2))
        cls.s = cls._genSignal(cls.size)[cls.size[0] // 2, :, :]
        for _ in range(100):
            cls.present.append(cls._genNoiseNPS(
                cls.size[:2], 2, 2) + cls.s * cls.contrast)

    def test_NPW(self):
        np.random.seed(7127433)
        MO = NPW(signal=self.s,
                 signal_present_samples=self.present,
                 signal_absent_samples=self.absent,
                 training_ratio=0.3)

        statistics, scores = MO.run_readers(10)

        self.assertTrue([s["auc"] for s in statistics] == self.truth["NPW"])

    def test_PWMF(self):
        np.random.seed(7127433)
        MO = PWMF(signal=self.s,
                  signal_present_samples=self.present,
                  signal_absent_samples=self.absent,
                  training_ratio=0.3)

        statistics, scores = MO.run_readers(10)

        self.assertTrue([s["auc"] for s in statistics] == self.truth["PWMF"])

    def test_CHO_LGauss(self):
        np.random.seed(7127433)
        MO = CHO(signal=self.s,
                 signal_present_samples=self.present,
                 signal_absent_samples=self.absent,
                 training_ratio=0.3,
                 channel_type="LGauss",
                 channel_params={
                              "LGauss_N": 5,
                              "LGauss_A": 25,
                              "LGauss_B": 25
                 })

        statistics, scores = MO.run_readers(10)

        self.assertTrue([s["auc"] for s in statistics] == self.truth["CHO_LG"])

    def test_CHO_DoG(self):
        np.random.seed(7127433)
        MO = CHO(signal=self.s,
                 signal_present_samples=self.present,
                 signal_absent_samples=self.absent,
                 training_ratio=0.3,
                 channel_type="DoG",
                 channel_params={
                     'DoG_N': 10,
                     'DoG_A': 1.67,
                     'DoG_Q': 1.4
                 })

        statistics, scores = MO.run_readers(10)

        self.assertTrue([s["auc"] for s in statistics] ==
                        self.truth["CHO_DoG"])

    def test_CHO_Gabor(self):
        np.random.seed(7127433)
        MO = CHO(signal=self.s,
                 signal_present_samples=self.present,
                 signal_absent_samples=self.absent,
                 training_ratio=0.3,
                 channel_type="Gabor",
                 channel_params={
                     'Gabor_PPD': 45,
                     'Gabor_O': 8,
                     'Gabor_F': 5
                 })

        statistics, scores = MO.run_readers(10)

        self.assertTrue([s["auc"] for s in statistics] ==
                        self.truth["CHO_Gabor"])

        pass

    # def test_isupper(self):
    #     self.assertTrue('FOO'.isupper())
    #     self.assertFalse('Foo'.isupper())

    # def test_split(self):
    #     s = 'hello world'
    #     self.assertEqual(s.split(), ['hello', 'world'])
    #     # check that s.split fails when the separator is not a string
    #     with self.assertRaises(TypeError):
    #         s.split(2)

    def _genNoiseNPS(size, noise_std, exponent):
        def _euclidean2D(dim):
            halfc = int(np.ceil((dim[1] - 1) / 2))
            halfr = int(np.ceil((dim[0] - 1) / 2))
            cols = np.linspace(halfc, -1 * (halfc - 1), dim[1])
            rows = np.linspace(halfr, -1 * (halfr - 1), dim[0])
            x, y = np.meshgrid(cols, rows)
            return [x, y]

        space = _euclidean2D(size)
        space[0] = space[0] / space[0].shape[0]
        space[1] = space[1] / space[1].shape[1]
        with np.errstate(divide='ignore'):
            NPS = np.fft.fftshift(
                1 / (np.sqrt(space[0] ** 2 + space[1]**2)**exponent))
        NPS[0, 0] = NPS[0, 1]
        white = np.random.rand(size[0], size[1]) * noise_std
        filtered = np.fft.ifftn(np.fft.fftn(white) * np.sqrt(NPS))
        return filtered.real / size[0]

    def _genSignal(size):
        def create_circular_mask3D(h, w, d, center=None, radius=None):
            if center is None:  # use the middle of the image
                center = (int(w / 2), int(h / 2), int(d / 2))
            if radius is None:  # use the smallest distance between the center and image walls
                radius = min(center[0], center[1], center[2], w -
                             center[0], h - center[1], d - center[2])

            Z, Y, X = np.ogrid[:d, :h, :w]
            dist_from_center = np.sqrt(
                (X - center[0])**2 + (Y - center[1])**2 + (Z - center[2])**2)

            mask = dist_from_center <= radius
            return mask
        s = create_circular_mask3D(
            size[0], size[1], size[2], radius=5).astype(np.float64)
        s = ndimage.gaussian_filter(s, 1)
        return s


if __name__ == '__main__':
    unittest.main()
    # t = TestModelObserver()
    # t.test_MO()
