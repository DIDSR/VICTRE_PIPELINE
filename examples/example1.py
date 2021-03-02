
#############################################################
#   MOVE THIS FILE TO THE PARENT DIRECTORY BEFORE RUNNING   #
#############################################################

from Victre import Pipeline
from Victre import Constants


seed = 1

roi_sizes = {Constants.VICTRE_SPICULATED: [65, 65, 5],
             Constants.VICTRE_CLUSTERCALC: [65, 65, 5]}

pline = Pipeline(seed=seed,
                 lesion_file="lesions/spiculated/mass_11_size1.00.h5",
                 roi_sizes=roi_sizes)

# we did not specify a pre-generated phantom, start from the beginning

pline.generate_phantom()

pline.compress_phantom(thickness=50)

pline.crop()

pline.insert_lesions(lesion_type=Constants.VICTRE_SPICULATED,
                     n=4)

pline.project()

pline.reconstruct()

pline.add_absent_ROIs(lesion_type=Constants.VICTRE_SPICULATED,
                      n=4)

pline.save_ROIs()
