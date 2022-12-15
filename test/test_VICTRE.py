import os
import sys
sys.path.append(os.path.abspath("."))  # noqa
from Victre import Pipeline, Constants
import unittest
import numpy as np
import hashlib


class TestPipeline(unittest.TestCase):
    roi_sizes = {Constants.VICTRE_SPICULATED: [65, 65, 5],
                 Constants.VICTRE_CLUSTERCALC: [65, 65, 5]}

    truth = {"generate": [505, 896, 717],
             "insert": {"phantom_coord": [[81, 76, 277, 2],
                                          [183, 29, 384, 2],
                                          [125, 138, 244, 2]],
                        "dbt_coord": [[616, 180, 16, 2],
                                      [364, 69, 37, 2],
                                      [694, 326, 25, 2]],
                        "dm_coord": [[191, 1480, 2],
                                     [76, 1204, 2],
                                     [350, 1563, 2]]
                        },
             "dbt_mask": "ebbc92604273bb53eab38abf1437579ee2462fdd686ed0b31f75384d6b096275"
             }

    def test_insert(self):
        # from new insertion
        plineI = Pipeline(seed=1,
                          phantom_file="phantoms/pc_1_crop.raw.gz",
                          lesion_file="lesions/spiculated/mass_11_size1.00.h5",
                          roi_sizes=self.roi_sizes,
                          results_folder="results_insert"
                          )
        np.random.seed(7127433)
        plineI.insert_lesions(lesion_type=Constants.VICTRE_SPICULATED,
                              n=3)
        for les in range(len(plineI.lesions)):
            for coord in range(len(plineI.lesions[les])):
                self.assertAlmostEqual(plineI.lesions[les][coord],
                                       np.array(self.truth["insert"]["phantom_coord"][les][coord]), delta=10)
            for coord in range(len(plineI.lesion_locations["dbt"][les])):
                self.assertAlmostEqual(plineI.lesion_locations["dbt"][les][coord],
                                       np.array(self.truth["insert"]["dbt_coord"][les][coord]), delta=10)
            for coord in range(len(plineI.lesion_locations["dm"][les])):
                self.assertAlmostEqual(plineI.lesion_locations["dm"][les][coord],
                                       np.array(self.truth["insert"]["dm_coord"][les][coord]), delta=10)

        # from loaded trial
        plineI = Pipeline(seed=1,
                          roi_sizes=self.roi_sizes,
                          results_folder="results_insert"
                          )

        for les in range(len(plineI.lesions)):
            for coord in range(len(plineI.lesions[les])):
                self.assertAlmostEqual(plineI.lesions[les][coord],
                                       np.array(self.truth["insert"]["phantom_coord"][les][coord]), delta=10)
            for coord in range(len(plineI.lesion_locations["dbt"][les])):
                self.assertAlmostEqual(plineI.lesion_locations["dbt"][les][coord],
                                       np.array(self.truth["insert"]["dbt_coord"][les][coord]), delta=10)
            for coord in range(len(plineI.lesion_locations["dm"][les])):
                self.assertAlmostEqual(plineI.lesion_locations["dm"][les][coord],
                                       np.array(self.truth["insert"]["dm_coord"][les][coord]), delta=10)

    def test_DBT_segmentation(self):
        plineD = Pipeline(seed=1,
                          phantom_file="phantoms/pc_1_crop.raw.gz",
                          lesion_file="lesions/spiculated/mass_11_size1.00.h5",
                          roi_sizes=self.roi_sizes,
                          results_folder="results_segmentation"
                          )
        segm = plineD.get_DBT_segmentation()
        self.assertEqual(hashlib.sha256(segm.tobytes()).hexdigest(),
                         self.truth["dbt_mask"])

    def test_generation(self):
        plineG = Pipeline(seed=7127433,
                          roi_sizes=self.roi_sizes,
                          results_folder="results_generation",
                          arguments_generation={
                              "imgRes": 0.1
                          }
                          )
        plineG.generate_phantom()
        self.assertAlmostEqual(
            plineG.arguments_mcgpu["number_voxels"], self.truth["generate"], delta=10)


if __name__ == '__main__':
    unittest.main(failfast=True)
    # t = TestPipeline()
    # t.setUpClass()
    # t.test_generation()
