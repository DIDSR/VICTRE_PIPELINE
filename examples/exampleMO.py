from Victre.ModelObserver import CHO, NPW, PWMF
import json
import numpy as np


results_folder = "results"  # change this to the VICTRE results folder
n_readers = 10  # number of random iterations for training and testing
training_ratio = 0.7  # by default testing ratio is 1 - training_ratio

# modalities to read from the results folder (dm and/or dbt)
modalities = ["dm", "dbt"]

# variables where results will be stored
statistics = {}
scores = {}


class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


for modality in modalities:
    # select the model we want to use

    MO = NPW(results_folder=results_folder,
             modality=modality,
             training_ratio=training_ratio)

    # MO = PWMF(results_folder=results_folder,
    #           modality=modality,
    #           training_ratio=training_ratio)

    # MO = CHO(channel_type="LGauss",
    #          channel_params={
    #              "LGauss_N": 5,
    #              "LGauss_A": 5,
    #              "LGauss_B": 5
    #          },
    #          results_folder=results_folder,
    #          modality=modality,
    #          training_ratio=training_ratio,
    #          channel_filter="convolutional")

    statistics[modality], scores[modality] = MO.run_readers(
        n_readers)

# saving results as iMRMC format
MO.save_iMRMC(results_folder + "output.imrmc",
              scores=scores)

# saving results as JSON format, including fpr, tpr, auc and d'
with open(results_folder + 'output.json', 'w') as f:
    json.dump(statistics, f, indent=4, cls=NumpyEncoder)

print(statistics)
