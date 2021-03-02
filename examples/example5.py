
#############################################################
#   MOVE THIS FILE TO THE PARENT DIRECTORY BEFORE RUNNING   #
#############################################################

from Victre.Pipeline import Pipeline
from Victre import Constants
import numpy as np


seed = 1

inserted_phantom = "phantoms/pcl_1.raw.gz"
roi_sizes = {Constants.VICTRE_SPICULATED: [65, 65, 5],
             Constants.VICTRE_CLUSTERCALC: [65, 65, 5]}

pline = Pipeline(seed=seed,
                 phantom_file=inserted_phantom,
                 roi_sizes=roi_sizes)

# phantom is ready to project, lesions were added in a previous step
pline.project()

pline.reconstruct()

pline.add_absent_ROIs(lesion_type=Constants.VICTRE_SPICULATED,
                      n=4)

pline.save_ROIs()
