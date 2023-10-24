import sys
from Victre.Pipeline import Pipeline
from Victre import Constants
import os
import argparse


parser = argparse.ArgumentParser(
    description='Runs on VICTRE pipeline task for DM vs DBT.')
parser.add_argument('--seed', type=int, default=1)
parser.add_argument('--density', type=str, default="dense")
parser.add_argument('--results', type=str,
                    help='Results folder.', default="./results")
parser.add_argument('--dose', type=float,
                    help='Number of histories on MCGPU.', default=2.04e9)
parser.add_argument('--simulate_lesions', action='store_true',
                    help='whether to simulate lesions')
parser.add_argument('--phantom_file', type=str, required=True)
parser.add_argument('--lesion_only', action='store_true')
parser.add_argument('--lesion_density', type=float,
                    help='lesion Density', default=1.1)
parser.add_argument('--lesion_file', type=str, default=None)
parser.add_argument('--flatfield_file', type=str, default="")

arguments = vars(parser.parse_args())

results_folder = arguments["results"]
seed = arguments["seed"]

if seed is None:
    seed = 1  # arguments["seed"]
density = arguments["density"]


lesion_file = arguments["lesion_file"]
flatfield_file = arguments["flatfield_file"]
singlePhantomFile = arguments["phantom_file"]

print('lesion file ' + str(lesion_file))

roi_sizes = {Constants.VICTRE_SPICULATED: [109, 109, 9],
             Constants.VICTRE_CLUSTERCALC: [65, 65, 5]}

histories = {"dense": 7.8e9,
             "hetero": 1.02e10,
             "scattered": 2.04e10,
             "fatty": 2.22e10}

LESION_SIMULATION = arguments["simulate_lesions"]
print("LESION_SIMULATION ", str(LESION_SIMULATION))

# else:
if arguments["dose"]:
    arguments_mcgpu = {"number_histories": arguments["dose"]}
else:
    arguments_mcgpu = {"number_histories": histories[arguments["density"]]}

print('density ', arguments["density"])
print('dose ', arguments_mcgpu["number_histories"])

spectrum_file = "./Victre/projection/spectrum/W28kVp_Rh50um_Be1mm.spc"


# INITIALIZE AND PROJECT
if "fatty" in arguments["density"]:
    spectrum_file = "./Victre/projection/spectrum/W30kVp_Rh50um_Be1mm.spc"
    arguments_mcgpu.update(
        {"fam_beam_aperture": [15.0, 11.2]})

# CONFIGURE PIPELINE FOR DM ONLY
arguments_mcgpu.update(
    {"number_projections": 1,
        "number_histories": arguments_mcgpu["number_histories"] * 25 * 2 / 3,
        "angular_blur": 0}
)

materials = Constants.VICTRE_DEFAULT_MATERIALS

if arguments["lesion_only"]:
    # set all other densities to 0, except lesion
    materials = [{"material": "./Victre/projection/material/air__5-120keV.mcgpu.gz",
                  "density": 0,
                  "voxel_id": list(Constants.PHANTOM_MATERIALS.values())
                  },
                 {"material": "./Victre/projection/material/glandular__5-120keV.mcgpu.gz",
                  "density": 1.06,
                  "voxel_id": [Constants.LESION_MATERIALS[Constants.VICTRE_SPICULATED]]
                  },
                 {"material": "./Victre/projection/material/CalciumOxalate__5-120keV.mcgpu.gz",
                  "density": 1.781,
                  "voxel_id": [Constants.LESION_MATERIALS[Constants.VICTRE_CLUSTERCALC]]
                  }
                 ]

os.makedirs(results_folder, exist_ok=True)

if os.path.isfile(flatfield_file):
    print('Using existing flatfield file!')
else:
    flatfield_file = None

if LESION_SIMULATION:
    # computing lesion density enhancement
    print('LESION DENSITY ' + str(arguments['lesion_density']))
    mass_mat = next((sub for sub in materials if sub["voxel_id"] == [
        Constants.LESION_MATERIALS[Constants.VICTRE_SPICULATED]]), None)
    mass_mat["density"] *= arguments['lesion_density']

    pline = Pipeline(results_folder=results_folder,
                     seed=int(seed),
                     phantom_file=singlePhantomFile,
                     spectrum_file=spectrum_file,
                     lesion_file=lesion_file,
                     flatfield_DM=flatfield_file,
                     roi_sizes=roi_sizes,
                     arguments_mcgpu=arguments_mcgpu,
                     materials=materials)

    print('inserting lesion!')
    # delete candidate locations for fully random insertion
    pline.candidate_locations = None
    pline.insert_lesions(
        lesion_type=Constants.VICTRE_SPICULATED, n=1)  # n=4


else:
    pline = Pipeline(results_folder=results_folder,
                     seed=int(seed),
                     phantom_file=singlePhantomFile,
                     spectrum_file=spectrum_file,
                     flatfield_DM=flatfield_file,
                     arguments_mcgpu=arguments_mcgpu,
                     materials=materials)

# perform projection
pline.project()
