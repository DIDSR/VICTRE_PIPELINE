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
        pline = Pipeline(seed=1,
                         phantom_file="phantoms/pc_1_crop.raw.gz",
                         lesion_file="lesions/spiculated/mass_11_size1.00.h5",
                         roi_sizes=self.roi_sizes,
                         results_folder="results_insert"
                         )
        np.random.seed(7127433)
        pline.insert_lesions(lesion_type=Constants.VICTRE_SPICULATED,
                             n=3)
        self.assertTrue(np.all(pline.lesions == np.array(
            self.truth["insert"]["phantom_coord"])))
        self.assertTrue(np.all(pline.lesion_locations["dbt"] == np.array(
            self.truth["insert"]["dbt_coord"])))
        self.assertTrue(np.all(pline.lesion_locations["dm"] == np.array(
            self.truth["insert"]["dm_coord"])))

    def test_DBT_segmentation(self):
        pline = Pipeline(seed=1,
                         phantom_file="phantoms/pc_1_crop.raw.gz",
                         lesion_file="lesions/spiculated/mass_11_size1.00.h5",
                         roi_sizes=self.roi_sizes,
                         results_folder="results_segmentation"
                         )
        segm = pline.get_DBT_segmentation()
        self.assertEqual(hashlib.sha256(segm.tobytes()).hexdigest(),
                         self.truth["dbt_mask"])

    def test_generation(self):
        pline = Pipeline(seed=7127433,
                         roi_sizes=self.roi_sizes,
                         results_folder="results_generation",
                         arguments_generation={
                             "imgRes": 0.1
                         }
                         )
        pline.generate_phantom()
        self.assertAlmostEqual(
            pline.arguments_mcgpu["number_voxels"], self.truth["generate"], delta=10)


if __name__ == '__main__':
    unittest.main()
    # t = TestPipeline()
    # t.setUpClass()
    # t.test_generation()
