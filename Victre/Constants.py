"""!
==================================================
        __   __ ___  ___  _____  ___  ___       
        \ \ / /|_ _|/ __||_   _|| _ \| __|      
         \ V /  | || (__   | |  |   /| _|       
          \_/  |___|\___|  |_|  |_|_\|___|      
     ___  ___  ___  ___  _     ___  _  _  ___ 
    | _ \|_ _|| _ \| __|| |   |_ _|| \| || __|
    |  _/ | | |  _/| _| | |__  | | | .` || _| 
    |_|  |___||_|  |___||____||___||_|\_||___|
    
==================================================

 Author: Miguel A. Lago
		 miguel.lago@fda.hhs.gov


				    DISCLAIMER

 This software and documentation (the "Software") were 
 developed at the Food and Drug Administration (FDA) by 
 employees of the Federal Government in the course of 
 their official duties. Pursuant to Title 17, Section 
 105 of the United States Code, this work is not subject
 to copyright protection and is in the public domain. 
 Permission is hereby granted, free of charge, to any 
 person obtaining a copy of the Software, to deal in the 
 Software without restriction, including without 
 limitation the rights to use, copy, modify, merge, 
 publish, distribute, sublicense, or sell copies of the 
 Software or derivatives, and to permit persons to whom 
 the Software is furnished to do so. FDA assumes no 
 responsibility whatsoever for use by other parties of 
 the Software, its source code, documentation or compiled 
 executables, and makes no guarantees, expressed or 
 implied, about its quality, reliability, or any other 
 characteristic. Further, use of this code in no way 
 implies endorsement by the FDA or confers any advantage 
 in regulatory decisions. Although this software can be 
 redistributed and/or modified freely, we ask that any 
 derivative works bear some notice that they are derived 
 from it, and any modified versions bear some notice that 
 they have been modified. 

 More information: https://github.com/DIDSR/VICTRE_PIPELINE

"""

VICTRE_CLUSTERCALC = 1
VICTRE_SPICULATED = 2

LESION_MATERIALS = {VICTRE_CLUSTERCALC: 250,
                    VICTRE_SPICULATED: 200}

PHANTOM_MATERIALS = {
    "air": 0,
    "adipose": 1,
    "skin": 2,
    "glandular": 29,
    "nipple": 33,
    "muscle": 40,
    "paddle": 50,
    "antiscatter_grid": 65,
    "detector": 66,
    "ligament": 88,
    "TDLU": 95,
    "duct": 125,
    "artery": 150,
    "vein": 225
}

INSERTION_MAX_TRIES = 2000
INSERTION_MAX_TOTAL_ATTEMPTS = 20

FLATFIELD_REPETITIONS = 1

FLATFIELD_DOSE_MULTIPLIER = 10

FORBIDDEN_OVERLAP = [PHANTOM_MATERIALS["air"], PHANTOM_MATERIALS["skin"], PHANTOM_MATERIALS["nipple"], PHANTOM_MATERIALS["muscle"]]

VICTRE_DEFAULT_MATERIALS = [
    {"material": "./Victre/projection/material/air__5-120keV.mcgpu.gz",
     "density": 0.0012,
     "voxel_id": [PHANTOM_MATERIALS["air"]]
     },
    {"material": "./Victre/projection/material/adipose__5-120keV.mcgpu.gz",
     "density": 0.92,
     "voxel_id": [PHANTOM_MATERIALS["adipose"]]
     },
    {"material": "./Victre/projection/material/skin__5-120keV.mcgpu.gz",
     "density": 1.09,
     "voxel_id": [PHANTOM_MATERIALS["skin"], PHANTOM_MATERIALS["nipple"]]
     },
    {"material": "./Victre/projection/material/glandular__5-120keV.mcgpu.gz",
     "density": 1.035,
     "voxel_id": [PHANTOM_MATERIALS["glandular"]]
     },
    {"material": "./Victre/projection/material/connective_Woodard__5-120keV.mcgpu.gz",
     "density": 1.12,
     "voxel_id": [PHANTOM_MATERIALS["ligament"]]
     },
    {"material": "./Victre/projection/material/muscle__5-120keV.mcgpu.gz",
     "density": 1.05,
     "voxel_id": [PHANTOM_MATERIALS["muscle"], PHANTOM_MATERIALS["duct"], PHANTOM_MATERIALS["TDLU"]]
     },
    {"material": "./Victre/projection/material/blood__5-120keV.mcgpu.gz",
     "density": 1.00,
     "voxel_id": [PHANTOM_MATERIALS["artery"], PHANTOM_MATERIALS["vein"]]
     },
    {"material": "./Victre/projection/material/polystyrene__5-120keV.mcgpu.gz",
     "density": 1.06,
     "voxel_id": [PHANTOM_MATERIALS["paddle"]]
     },
    {"material": "./Victre/projection/material/glandular__5-120keV.mcgpu.gz",  # MASS MATERIAL
     # "density": 1.06, # VICTRE 1.0
     "density": 1.06,
     "voxel_id": [LESION_MATERIALS[VICTRE_SPICULATED]]
     },
    {"material": "./Victre/projection/material/CalciumOxalate__5-120keV.mcgpu.gz",  # CLUSTER CALC MATERIAL
     "density": 1.781,
     "voxel_id": [LESION_MATERIALS[VICTRE_CLUSTERCALC]]
     },
    {"material": "./Victre/projection/material/W__5-120keV.mcgpu.gz",
     "density": 1.0,
     "voxel_id": [PHANTOM_MATERIALS["antiscatter_grid"]]
     },
    {"material": "./Victre/projection/material/Se__5-120keV.mcgpu.gz",
     "density": 1.0,
     "voxel_id": [PHANTOM_MATERIALS["detector"]]
     },
]

VICTRE_DEFAULT_MCGPU = {
    "number_histories": 7.8e9,
    "random_seed": 31415990,
    "selected_gpu": 0,
    "number_gpus": 1,
    "gpu_threads": 128,
    "histories_per_thread": 5000,
    "spectrum_file": "./Victre/projection/spectrum/W28kVp_Rh50um_Be1mm.spc",
    "source_position": [0.00001, 4.825, 63.0],
    "source_direction": [0.0, 0.0, -1.0],
    "fam_beam_aperture": [15.0, 7.4686667],
    "euler_angles": [90.0, -90.0, 180.0],
    "focal_spot": 0.0300,
    "angular_blur": 0.18,
    "collimate_beam": "YES",
    "output_file": None,
    "image_pixels": [3000, 1500],
    "image_size": [25.50, 12.75],
    "distance_source": 65.00,
    "image_offset": [0, 0],
    "detector_thickness": 0.02,
    "mean_free_path": 0.004027,
    "k_edge_energy": [12658.0, 11223.0, 0.596, 0.00593],
    "detector_gain": [50.0, 0.99],
    "additive_noise": 5200.0,
    "cover_thickness": [0.10, 1.9616],
    "antiscatter_grid_ratio": [5.0, 31.0, 0.0065],
    "antiscatter_strips": [0.00089945, 1.9616],
    "antiscatter_grid_lines": 0,
    "number_projections": 25,
    "rotation_axis_distance": 60.0,
    "projections_angle": 2.083333333333,
    "angular_rotation_first": -25.0,
    "rotation_axis": [1.0, 0.0, 0.0],
    "axis_translation": 0,
    "detector_fixed": "YES",
    "simulate_both": "YES",
    "tally_material_dose": "YES",
    "tally_voxel_dose": "NO",
    "output_dose_filename": "mc-gpu_dose.dat",
    "roi_voxel_dose_x": [1, 751],
    "roi_voxel_dose_y": [1, 1301],
    "roi_voxel_dose_z": [250, 250],
    "phantom_file": None,
    "voxel_geometry_offset": [0, 0, 0],
    "number_voxels": [810, 1920, 745],
    "voxel_size": [0.005, 0.005, 0.005],
    "low_resolution_voxel_size": [0, 0, 0]
}

VICTRE_DENSE = {
    "compressionThickness": 35,  # mm
    # phantom voxel size (mm)
    "imgRes": 0.05,
    # thickness of breast skin (mm)
    "skinThick": 0.75,
    # nipple length (mm)
    "nippleLen": 4.0,
    # nipple radius (mm)
    "nippleRad": 4.0,
    # nipple radius (mm)
    "areolaRad": 8.0,
    # left breast - select left or right breast (boolean)
    "leftBreast": "false",
    # desired fat fraction
    "targetFatFrac": 0.4,

    #####################
    # breast surface shape
    #####################
    # u resolution of base shape
    "surface_ures": 0.005,
    # v resolution of base shape
    "surface_vres": 0.005,
    # minimum point separation (mm)
    "surface_pointSep": 0.005,
    # back ring thickness (mm)
    "surface_ringWidth": 10.0,
    # back ring step size (mm)
    "surface_ringSep": 0.5,
    # angle to preserve while smoothing (degrees)
    "surface_featureAngle": 20.0,
    # fraction of triangles to decimate
    "surface_targetReduction": 0.05,
    # bottom scale was 1.3,
    "surface_a1b": 1.0,
    # top scale was 1.3,
    "surface_a1t": 1.0,
    # left scale
    "surface_a2l": 1.0,
    # right scale
    "surface_a2r": 1.0,
    # outward scale was 1.75,
    "surface_a3": 1.0,
    # u quadric exponent
    "surface_eps1": 1.2,
    # v quadric exponent
    "surface_eps2": 1.0,
    # do ptosis deformation (boolean)
    "surface_doPtosis": "true",
    "surface_ptosisB0": 0.2,
    "surface_ptosisB1": 0.05,
    # do turn deformation (boolean)
    "surface_doTurn": "false",
    "surface_turnC0": -0.498,
    "surface_turnC1": 0.213,
    # do top shape deformation (boolean)
    "surface_doTopShape": "true",
    "surface_topShapeS0": 0.0,
    "surface_topShapeS1": 0.0,
    "surface_topShapeT0": -12.0,
    "surface_topShapeT1": -5.0,
    # do flatten size deformation (boolean)
    "surface_doFlattenSide": "true",
    "surface_flattenSideG0": 1.5,
    "surface_flattenSideG1": -0.5,
    # do turn top deformation (boolean)
    "surface_doTurnTop": "true",
    "surface_turnTopH0": 0.166,
    "surface_turnTopH1": -0.372,

    #####################
    # breast compartment
    #####################
    # number of breast compartments
    "compartment_num": 10,
    # distance along nipple line of compartment seed base (mm)
    "compartment_seedBaseDist": 16,
    # fraction of phantom in nipple direction forced to be fat
    "compartment_backFatBufferFrac": 0.008,
    # number of backplane seed points
    "compartment_numBackSeeds": 150,
    # maximum seed jitter (fraction of subtended angle)
    "compartment_angularJitter": 0.125,
    # maximum seed jitter in nipple direction (mm)
    "compartment_zJitter": 5.0,
    # maximum radial distance from base seed as a fraction of distance to breast surface
    "compartment_maxFracRadialDist": 0.5,
    # minimum radial distance from base seed as a fraction of distance to breast surface
    "compartment_minFracRadialDist": 0.25,
    # minimum scale in nipple direction
    "compartment_minScaleNippleDir": 0.01,
    # maximum scale in nipple direction
    "compartment_maxScaleNippleDir": 0.01,
    # minimum scale in non-nipple direction
    "compartment_minScale": 30.0,
    # maximum scale in non-nipple direction
    "compartment_maxScale": 40.0,
    # minimum gland strength
    "compartment_minGlandStrength": 30.0,
    # maximum gland strength
    "compartment_maxGlandStrength": 30.0,
    # maximum compartment deflection angle from pointing towards nipple (fraction of pi)
    "compartment_maxDeflect": 0.01,
    # minimum scale skin seeds in nipple direction
    "compartment_minSkinScaleNippleDir": 10.0,
    # maximum scale skin seeds in nipple direction
    "compartment_maxSkinScaleNippleDir": 20.0,
    # minimum scale skin in non-nipple direction
    "compartment_minSkinScale": 200.0,
    # maximum scale skin in non-nipple direction
    "compartment_maxSkinScale": 400.0,
    # skin strength
    "compartment_skinStrength": 0.5,
    # back scale
    "compartment_backScale": 60.0,
    # back strength
    "compartment_backStrength": 1.0,
    # nipple scale
    "compartment_nippleScale": 5.0,
    # nipple strength
    "compartment_nippleStrength": 10.0,
    # check seeds within radius (mm)
    "compartment_voronSeedRadius": 100.0,

    #####################
    # TDLU variables
    #####################
    # maximum TDLU length
    "tdlu_maxLength": 2.0,
    # minimum TDLU length
    "tdlu_minLength": 1.0,
    # maximum TDLU width
    "tdlu_maxWidth": 1.0,
    # minimum TDLU width
    "tdlu_minWidth": 0.5,

    #####################
    # Perlin noise variables
    #####################
    # maximum fraction of radius deviation
    "perlin_noise_maxDeviation": 0.1,
    # starting frequency
    "perlin_noise_frequency": 0.1,
    # octave frequency multiplier
    "perlin_noise_lacunarity": 2.0,
    # octave signal decay
    "perlin_noise_persistence": 0.5,
    # number of frequency octaves
    "perlin_noise_numOctaves": 6,
    # x direction noise generation seed
    "perlin_noise_xNoiseGen": 683,
    # y direction noise generation seed
    "perlin_noise_yNoiseGen": 4933,
    # z direction noise generation seed
    "perlin_noise_zNoiseGen": 23,
    # seed noise generation
    "perlin_noise_seedNoiseGen": 3095,
    # shift noise generation seed
    "perlin_noise_shiftNoiseGen": 11,

    #####################
    # Compartment boundary noise
    #####################
    # maximum fraction of distance deviation
    "compartment_boundary_maxDeviation": 0.1,
    # starting frequency
    "compartment_boundary_frequency": 0.15,
    # octave frequency multiplier
    "compartment_boundary_lacunarity": 1.5,
    # octave signal decay
    "compartment_boundary_persistence": 0.5,

    #####################
    # Lobule boundary perturbation noise
    #####################
    # maximum fraction of distance deviation
    "lobule_boundary_maxDeviation": 0.25,
    # starting frequency
    "lobule_boundary_frequency": 0.09,
    # octave frequency multiplier
    "lobule_boundary_lacunarity": 2.0,
    # octave signal decay
    "lobule_boundary_persistence": 0.4,

    #####################
    # Lobule glandular buffer noise
    #####################
    # maximum fraction of distance deviation
    "lobule_glandular_maxDeviation": 0.15,
    # starting frequency
    "lobule_glandular_frequency": 0.05,
    # octave frequency multiplier
    "lobule_glandular_lacunarity": 1.5,
    # octave signal decay
    "lobule_glandular_persistence": 0.5,

    #####################
    # Voronoi segmentation variables
    #####################
    # fat voronoi seed density (mm^-3,)
    "voronoi_fatInFatSeedDensity": 0.001,
    # fat voronoi seed in glandular tissue density (mm^-3,)
    "voronoi_fatInGlandSeedDensity": 0.001,
    # glandular voronoi seed density (mm^-3,)
    "voronoi_glandInGlandSeedDensity": 0.0005,
    # maximum deflection (fraction of pi)
    "voronoi_TDLUDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenTDLU": 0.1,
    # maximum length scale
    "voronoi_maxScaleLenTDLU": 0.2,
    # minimum width scale
    "voronoi_minScaleWidTDLU": 40.0,
    # maximum width scale
    "voronoi_maxScaleWidTDLU": 45.0,
    # minimum strength
    "voronoi_minStrTDLU": 20.0,
    # maximum strength
    "voronoi_maxStrTDLU": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInFatDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInFat": 5.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInFat": 10.0,
    # minimum width scale
    "voronoi_minScaleWidFatInFat": 50.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInFat": 60.0,
    # minimum strength
    "voronoi_minStrFatInFat": 40.0,
    # maximum strength
    "voronoi_maxStrFatInFat": 50.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidFatInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInGland": 40.0,
    # minimum strength
    "voronoi_minStrFatInGland": 20.0,
    # maximum strength
    "voronoi_maxStrFatInGland": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_glandInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenGlandInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenGlandInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidGlandInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidGlandInGland": 40.0,
    # minimum strength
    "voronoi_minStrGlandInGland": 20.0,
    # maximum strength
    "voronoi_maxStrGlandInGland": 22.0,
    # check seeds in radius (mm)
    "voronoi_seedRadius": 40.0,

    #####################
    # fat variables
    #####################
    # min lobule axis length (mm)
    # was 20,
    "fat_minLobuleAxis": 20.0,
    # max lobule axis length (mm)
    # was 35,
    "fat_maxLobuleAxis": 30.0,
    # axial ratio min
    "fat_minAxialRatio": 0.13,
    # axial ratio max
    "fat_maxAxialRatio": 0.75,
    # minimum ligament separation between lobules
    "fat_minLobuleGap": 0.15,
    # maximum of absolute value of Fourier coefficient as fraction of main radius
    "fat_maxCoeffStr": 0.1,
    # minimum of absolute value of Fourier coefficient as fraction of main radius
    "fat_minCoeffStr": 0.05,
    # maximum number of trial lobules
    # was 401,
    "fat_maxLobuleTry": 401,


    #####################
    # ligament variables
    #####################
    "ligament_thickness": 0.1,
    "ligament_targetFrac": 0.85,
    "ligament_maxTry": 15000,
    "ligament_minAxis": 20.0,
    "ligament_maxAxis": 25.0,
    "ligament_minAxialRatio": 0.2,
    "ligament_maxAxialRatio": 0.3,
    "ligament_maxPerturb": 0.05,
    "ligament_maxDeflect": 0.12,
    "ligament_scale": 0.007,
    "ligament_lacunarity": 1.5,
    "ligament_persistence": 0.3,
    "ligament_numOctaves": 6,

    #####################
    # duct tree variables
    #####################
    # target number of branches (uint)
    "duct_tree_maxBranch": 400,
    # maximum generation (uint)
    "duct_tree_maxGen": 7,
    # initial radius of tree (mm)
    "duct_tree_initRad": 0.5,
    # base Length of root duct at nipple (mm)
    "duct_tree_baseLength": 7.6,
    # number of voxels for tree density tracking (uint)
    "duct_tree_nFillX": 50,
    "duct_tree_nFillY": 50,
    "duct_tree_nFillZ": 50,

    #####################
    # duct branch variables
    #####################
    # minimum branch radius to have children (mm)
    "duct_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "duct_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "duct_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "duct_branch_lenShrink": 0.5,
    # maximum jitter in branch length (fraction)
    "duct_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "duct_branch_rotateJitter": 0.1,

    #####################
    # duct segment variables
    #####################
    # radius distribution shape parameters
    "duct_segment_radiusBetaA": 6.0,
    "duct_segment_radiusBetaB": 10.0,
    # fraction of branch length per segment
    "duct_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "duct_segment_maxCurvRad": 10.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "duct_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "duct_segment_minEndRad": 0.95,
    "duct_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "duct_segment_angleWt": 1.0,
    # cost function density weighting
    "duct_segment_densityWt": 20.0,
    # number of trial segments to generate (uint)
    "duct_segment_numTry": 50,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "duct_segment_maxTry": 100,
    # total number of segment tries before completely giving up
    "duct_segment_absMaxTry": 10000,
    # step size for checking segment is valid (mm)
    "duct_segment_roiStep": 0.1,


    #####################
    # vessel tree variables
    #####################
    # target number of branches (uint)
    "vessel_tree_maxBranch": 750,
    # maximum generation (uint)
    "vessel_tree_maxGen": 6,
    # initial radius of tree (mm)
    "vessel_tree_initRad": 0.75,
    # base length of root vessel (mm)
    "vessel_tree_baseLength": 15.0,
    # number of voxels for tree density tracking (uint)
    "vessel_tree_nFillX": 30,
    "vessel_tree_nFillY": 69,
    "vessel_tree_nFillZ": 69,

    #####################
    # vessel branch variables
    #####################
    # minimum branch radius to have children (mm)
    "vessel_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "vessel_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "vessel_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "vessel_branch_lenShrink": 0.8,
    # maximum jitter in branch length (fraction)
    "vessel_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "vessel_branch_rotateJitter": 0.1,

    #####################
    # vessel segment variables
    #####################
    # radius distribution shape parameters
    "vessel_segment_radiusBetaA": 6.0,
    "vessel_segment_radiusBetaB": 10.0,
    # fraction of branch length to segment
    "vessel_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "vessel_segment_maxCurvRad": 200.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "vessel_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "vessel_segment_minEndRad": 0.95,
    "vessel_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "vessel_segment_angleWt": 100.0,
    # cost function density weighting
    "vessel_segment_densityWt": 1.0,
    # cost function direction weighting
    "vessel_segment_dirWt": 100.0,
    # number of trial segments to generate (uint)
    "vessel_segment_numTry": 100,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "vessel_segment_maxTry": 300,
    # total number of segment tries before completely giving up
    "vessel_segment_absMaxTry": 100000,
    # step size for checking segment is valid (mm)
    "vessel_segment_roiStep": 0.1,
}

VICTRE_HETERO = {
    "compressionThickness": 45,  # mm
    # phantom voxel size (mm)
    "imgRes": 0.05,
    # thickness of breast skin (mm)
    "skinThick": 0.75,
    # nipple length (mm)
    "nippleLen": 4.0,
    # nipple radius (mm)
    "nippleRad": 4.0,
    # nipple radius (mm)
    "areolaRad": 8.0,
    # left breast - select left or right breast (boolean)
    "leftBreast": "false",
    # desired fat fraction
    "targetFatFrac": 0.66,

    #####################
    # breast surface shape
    #####################
    # u resolution of base shape
    "surface_ures": 0.005,
    # v resolution of base shape
    "surface_vres": 0.005,
    # minimum point separation (mm)
    "surface_pointSep": 0.005,
    # back ring thickness (mm)
    "surface_ringWidth": 10.0,
    # back ring step size (mm)
    "surface_ringSep": 0.5,
    # angle to preserve while smoothing (degrees)
    "surface_featureAngle": 20.0,
    # fraction of triangles to decimate
    "surface_targetReduction": 0.05,
    # bottom scale was 1.3,
    "surface_a1b": 1.2,
    # top scale was 1.3,
    "surface_a1t": 1.2,
    # left scale
    "surface_a2l": 1.0,
    # right scale
    "surface_a2r": 1.0,
    # outward scale was 1.75,
    "surface_a3": 1.55,
    # u quadric exponent
    "surface_eps1": 1.2,
    # v quadric exponent
    "surface_eps2": 1.0,
    # do ptosis deformation (boolean)
    "surface_doPtosis": "true",
    "surface_ptosisB0": 0.2,
    "surface_ptosisB1": 0.05,
    # do turn deformation (boolean)
    "surface_doTurn": "false",
    "surface_turnC0": -0.498,
    "surface_turnC1": 0.213,
    # do top shape deformation (boolean)
    "surface_doTopShape": "true",
    "surface_topShapeS0": 0.0,
    "surface_topShapeS1": 0.0,
    "surface_topShapeT0": -12.0,
    "surface_topShapeT1": -5.0,
    # do flatten size deformation (boolean)
    "surface_doFlattenSide": "true",
    "surface_flattenSideG0": 1.5,
    "surface_flattenSideG1": -0.5,
    # do turn top deformation (boolean)
    "surface_doTurnTop": "true",
    "surface_turnTopH0": 0.166,
    "surface_turnTopH1": -0.372,

    #####################
    # breast compartment
    #####################
    # number of breast compartments
    "compartment_num": 10,
    # distance along nipple line of compartment seed base (mm)
    "compartment_seedBaseDist": 16,
    # fraction of phantom in nipple direction forced to be fat
    "compartment_backFatBufferFrac": 0.01,
    # number of backplane seed points
    "compartment_numBackSeeds": 250,
    # maximum seed jitter (fraction of subtended angle)
    "compartment_angularJitter": 0.125,
    # maximum seed jitter in nipple direction (mm)
    "compartment_zJitter": 5.0,
    # maximum radial distance from base seed as a fraction of distance to breast surface
    "compartment_maxFracRadialDist": 0.5,
    # minimum radial distance from base seed as a fraction of distance to breast surface
    "compartment_minFracRadialDist": 0.25,
    # minimum scale in nipple direction
    "compartment_minScaleNippleDir": 0.01,
    # maximum scale in nipple direction
    "compartment_maxScaleNippleDir": 0.01,
    # minimum scale in non-nipple direction
    "compartment_minScale": 30.0,
    # maximum scale in non-nipple direction
    "compartment_maxScale": 40.0,
    # minimum gland strength
    "compartment_minGlandStrength": 30.0,
    # maximum gland strength
    "compartment_maxGlandStrength": 30.0,
    # maximum compartment deflection angle from pointing towards nipple (fraction of pi)
    "compartment_maxDeflect": 0.01,
    # minimum scale skin seeds in nipple direction
    "compartment_minSkinScaleNippleDir": 10.0,
    # maximum scale skin seeds in nipple direction
    "compartment_maxSkinScaleNippleDir": 20.0,
    # minimum scale skin in non-nipple direction
    "compartment_minSkinScale": 200.0,
    # maximum scale skin in non-nipple direction
    "compartment_maxSkinScale": 200.0,
    # skin strength
    "compartment_skinStrength": 1.0,
    # back scale
    "compartment_backScale": 60.0,
    # back strength
    "compartment_backStrength": 2.0,
    # nipple scale
    "compartment_nippleScale": 5.0,
    # nipple strength
    "compartment_nippleStrength": 10.0,
    # check seeds within radius (mm)
    "compartment_voronSeedRadius": 100.0,

    #####################
    # TDLU variables
    #####################
    # maximum TDLU length
    "tdlu_maxLength": 2.0,
    # minimum TDLU length
    "tdlu_minLength": 1.0,
    # maximum TDLU width
    "tdlu_maxWidth": 1.0,
    # minimum TDLU width
    "tdlu_minWidth": 0.5,

    #####################
    # Perlin noise variables
    #####################
    # maximum fraction of radius deviation
    "perlin_noise_maxDeviation_perlin": 0.1,
    # starting frequency
    "perlin_noise_frequency_perlin": 0.1,
    # octave frequency multiplier
    "perlin_noise_lacunarity_perlin": 2.0,
    # octave signal decay
    "perlin_noise_persistence_perlin": 0.5,
    # number of frequency octaves
    "perlin_noise_numOctaves": 6,
    # x direction noise generation seed
    "perlin_noise_xNoiseGen": 683,
    # y direction noise generation seed
    "perlin_noise_yNoiseGen": 4933,
    # z direction noise generation seed
    "perlin_noise_zNoiseGen": 23,
    # seed noise generation
    "perlin_noise_seedNoiseGen": 3095,
    # shift noise generation seed
    "perlin_noise_shiftNoiseGen": 11,

    #####################
    # Compartment boundary noise
    #####################
    # maximum fraction of distance deviation
    "compartment_boundary_maxDeviation": 0.1,
    # starting frequency
    "compartment_boundary_frequency": 0.15,
    # octave frequency multiplier
    "compartment_boundary_lacunarity": 1.5,
    # octave signal decay
    "compartment_boundary_persistence": 0.5,

    #####################
    # Lobule boundary perturbation noise
    #####################
    # maximum fraction of distance deviation
    "lobule_boundary_maxDeviation": 0.25,
    # starting frequency
    "lobule_boundary_frequency": 0.09,
    # octave frequency multiplier
    "lobule_boundary_lacunarity": 2.0,
    # octave signal decay
    "lobule_boundary_persistence": 0.4,

    #####################
    # Lobule glandular buffer noise
    #####################
    # maximum fraction of distance deviation
    "lobule_glandular_maxDeviation": 0.15,
    # starting frequency
    "lobule_glandular_frequency": 0.05,
    # octave frequency multiplier
    "lobule_glandular_lacunarity": 1.5,
    # octave signal decay
    "lobule_glandular_persistence": 0.5,

    #####################
    # Voronoi segmentation variables
    #####################
    # fat voronoi seed density (mm^-3,)
    "voronoi_fatInFatSeedDensity": 0.001,
    # fat voronoi seed in glandular tissue density (mm^-3,)
    "voronoi_fatInGlandSeedDensity": 0.001,
    # glandular voronoi seed density (mm^-3,)
    "voronoi_glandInGlandSeedDensity": 0.0005,
    # maximum deflection (fraction of pi)
    "voronoi_TDLUDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenTDLU": 0.1,
    # maximum length scale
    "voronoi_maxScaleLenTDLU": 0.2,
    # minimum width scale
    "voronoi_minScaleWidTDLU": 40.0,
    # maximum width scale
    "voronoi_maxScaleWidTDLU": 45.0,
    # minimum strength
    "voronoi_minStrTDLU": 20.0,
    # maximum strength
    "voronoi_maxStrTDLU": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInFatDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInFat": 5.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInFat": 10.0,
    # minimum width scale
    "voronoi_minScaleWidFatInFat": 50.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInFat": 60.0,
    # minimum strength
    "voronoi_minStrFatInFat": 40.0,
    # maximum strength
    "voronoi_maxStrFatInFat": 50.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidFatInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInGland": 40.0,
    # minimum strength
    "voronoi_minStrFatInGland": 20.0,
    # maximum strength
    "voronoi_maxStrFatInGland": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_glandInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenGlandInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenGlandInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidGlandInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidGlandInGland": 40.0,
    # minimum strength
    "voronoi_minStrGlandInGland": 20.0,
    # maximum strength
    "voronoi_maxStrGlandInGland": 22.0,
    # check seeds in radius (mm)
    "voronoi_seedRadius": 40.0,

    #####################
    # fat variables
    #####################
    # min lobule axis length (mm)
    # was 20,
    "fat_minLobuleAxis": 20.0,
    # max lobule axis length (mm)
    # was 35,
    "fat_maxLobuleAxis": 30.0,
    # axial ratio min
    "fat_minAxialRatio": 0.13,
    # axial ratio max
    "fat_maxAxialRatio": 0.75,
    # minimum ligament separation between lobules
    "fat_minLobuleGap": 0.15,
    # maximum of absolute value of Fourier coefficient as fraction of main radius
    "fat_maxCoeffStr": 0.1,
    # minimum of absolute value of Fourier coefficient as fraction of main radius
    "fat_minCoeffStr": 0.05,
    # maximum number of trial lobules
    # was 401,
    "fat_maxLobuleTry": 401,


    #####################
    # ligament variables
    #####################
    "ligament_thickness": 0.1,
    "ligament_targetFrac": 0.85,
    "ligament_maxTry": 15000,
    "ligament_minAxis": 20.0,
    "ligament_maxAxis": 25.0,
    "ligament_minAxialRatio": 0.2,
    "ligament_maxAxialRatio": 0.3,
    "ligament_maxPerturb": 0.05,
    "ligament_maxDeflect": 0.12,
    "ligament_scale": 0.007,
    "ligament_lacunarity": 1.5,
    "ligament_persistence": 0.3,
    "ligament_numOctaves": 6,

    #####################
    # duct tree variables
    #####################
    # target number of branches (uint)
    "duct_tree_maxBranch": 400,
    # maximum generation (uint)
    "duct_tree_maxGen": 7,
    # initial radius of tree (mm)
    "duct_tree_initRad": 0.5,
    # base Length of root duct at nipple (mm)
    "duct_tree_baseLength": 7.6,
    # number of voxels for tree density tracking (uint)
    "duct_tree_nFillX": 50,
    "duct_tree_nFillY": 50,
    "duct_tree_nFillZ": 50,

    #####################
    # duct branch variables
    #####################
    # minimum branch radius to have children (mm)
    "duct_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "duct_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "duct_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "duct_branch_lenShrink": 0.5,
    # maximum jitter in branch length (fraction)
    "duct_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "duct_branch_rotateJitter": 0.1,

    #####################
    # duct segment variables
    #####################
    # radius distribution shape parameters
    "duct_segment_radiusBetaA": 6.0,
    "duct_segment_radiusBetaB": 10.0,
    # fraction of branch length per segment
    "duct_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "duct_segment_maxCurvRad": 10.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "duct_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "duct_segment_minEndRad": 0.95,
    "duct_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "duct_segment_angleWt": 1.0,
    # cost function density weighting
    "duct_segment_densityWt": 20.0,
    # number of trial segments to generate (uint)
    "duct_segment_numTry": 50,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "duct_segment_maxTry": 100,
    # total number of segment tries before completely giving up
    "duct_segment_absMaxTry": 10000,
    # step size for checking segment is valid (mm)
    "duct_segment_roiStep": 0.1,


    #####################
    # vessel tree variables
    #####################
    # target number of branches (uint)
    "vessel_tree_maxBranch": 750,
    # maximum generation (uint)
    "vessel_tree_maxGen": 6,
    # initial radius of tree (mm)
    "vessel_tree_initRad": 0.75,
    # base length of root vessel (mm)
    "vessel_tree_baseLength": 15.0,
    # number of voxels for tree density tracking (uint)
    "vessel_tree_nFillX": 30,
    "vessel_tree_nFillY": 69,
    "vessel_tree_nFillZ": 69,

    #####################
    # vessel branch variables
    #####################
    # minimum branch radius to have children (mm)
    "vessel_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "vessel_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "vessel_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "vessel_branch_lenShrink": 0.8,
    # maximum jitter in branch length (fraction)
    "vessel_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "vessel_branch_rotateJitter": 0.1,

    #####################
    # vessel segment variables
    #####################
    # radius distribution shape parameters
    "vessel_segment_radiusBetaA": 6.0,
    "vessel_segment_radiusBetaB": 10.0,
    # fraction of branch length to segment
    "vessel_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "vessel_segment_maxCurvRad": 200.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "vessel_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "vessel_segment_minEndRad": 0.95,
    "vessel_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "vessel_segment_angleWt": 100.0,
    # cost function density weighting
    "vessel_segment_densityWt": 1.0,
    # cost function direction weighting
    "vessel_segment_dirWt": 100.0,
    # number of trial segments to generate (uint)
    "vessel_segment_numTry": 100,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "vessel_segment_maxTry": 300,
    # total number of segment tries before completely giving up
    "vessel_segment_absMaxTry": 100000,
    # step size for checking segment is valid (mm)
    "vessel_segment_roiStep": 0.1,
}

VICTRE_SCATTERED = {
    "compressionThickness": 55,  # mm
    # phantom voxel size (mm)
    "imgRes": 0.05,
    # thickness of breast skin (mm)
    "skinThick": 0.75,
    # nipple length (mm)
    "nippleLen": 4.0,
    # nipple radius (mm)
    "nippleRad": 4.0,
    # nipple radius (mm)
    "areolaRad": 8.0,
    # left breast - select left or right breast (boolean)
    "leftBreast": "false",
    # desired fat fraction
    "targetFatFrac": 0.85,

    #####################
    # breast surface shape
    #####################
    # u resolution of base shape
    "surface_ures": 0.005,
    # v resolution of base shape
    "surface_vres": 0.005,
    # minimum point separation (mm)
    "surface_pointSep": 0.005,
    # back ring thickness (mm)
    "surface_ringWidth": 10.0,
    # back ring step size (mm)
    "surface_ringSep": 0.5,
    # angle to preserve while smoothing (degrees)
    "surface_featureAngle": 20.0,
    # fraction of triangles to decimate
    "surface_targetReduction": 0.05,
    # bottom scale was 1.3,
    "surface_a1b": 1.4,
    # top scale was 1.3,
    "surface_a1t": 1.4,
    # left scale
    "surface_a2l": 1.2,
    # right scale
    "surface_a2r": 1.2,
    # outward scale was 1.75,
    "surface_a3": 2.0,
    # u quadric exponent
    "surface_eps1": 1.2,
    # v quadric exponent
    "surface_eps2": 1.0,
    # do ptosis deformation (boolean)
    "surface_doPtosis": "true",
    "surface_ptosisB0": 0.2,
    "surface_ptosisB1": 0.05,
    # do turn deformation (boolean)
    "surface_doTurn": "false",
    "surface_turnC0": -0.498,
    "surface_turnC1": 0.213,
    # do top shape deformation (boolean)
    "surface_doTopShape": "true",
    "surface_topShapeS0": 0.0,
    "surface_topShapeS1": 0.0,
    "surface_topShapeT0": -8.0,
    "surface_topShapeT1": -2.0,
    # do flatten size deformation (boolean)
    "surface_doFlattenSide": "true",
    "surface_flattenSideG0": 1.5,
    "surface_flattenSideG1": -0.5,
    # do turn top deformation (boolean)
    "surface_doTurnTop": "true",
    "surface_turnTopH0": 0.166,
    "surface_turnTopH1": -0.372,

    #####################
    # breast compartment
    #####################
    # number of breast compartments
    "compartment_num": 10,
    # distance along nipple line of compartment seed base (mm)
    "compartment_seedBaseDist": 16,
    # fraction of phantom in nipple direction forced to be fat
    "compartment_backFatBufferFrac": 0.01,
    # number of backplane seed points
    "compartment_numBackSeeds": 250,
    # maximum seed jitter (fraction of subtended angle)
    "compartment_angularJitter": 0.125,
    # maximum seed jitter in nipple direction (mm)
    "compartment_zJitter": 5.0,
    # maximum radial distance from base seed as a fraction of distance to breast surface
    "compartment_maxFracRadialDist": 0.5,
    # minimum radial distance from base seed as a fraction of distance to breast surface
    "compartment_minFracRadialDist": 0.25,
    # minimum scale in nipple direction
    "compartment_minScaleNippleDir": 0.01,
    # maximum scale in nipple direction
    "compartment_maxScaleNippleDir": 0.01,
    # minimum scale in non-nipple direction
    "compartment_minScale": 30.0,
    # maximum scale in non-nipple direction
    "compartment_maxScale": 40.0,
    # minimum gland strength
    "compartment_minGlandStrength": 30.0,
    # maximum gland strength
    "compartment_maxGlandStrength": 30.0,
    # maximum compartment deflection angle from pointing towards nipple (fraction of pi)
    "compartment_maxDeflect": 0.01,
    # minimum scale skin seeds in nipple direction
    "compartment_minSkinScaleNippleDir": 5.0,
    # maximum scale skin seeds in nipple direction
    "compartment_maxSkinScaleNippleDir": 5.0,
    # minimum scale skin in non-nipple direction
    "compartment_minSkinScale": 200.0,
    # maximum scale skin in non-nipple direction
    "compartment_maxSkinScale": 200.0,
    # skin strength
    "compartment_skinStrength": 2.0,
    # back scale
    "compartment_backScale": 60.0,
    # back strength
    "compartment_backStrength": 4.0,
    # nipple scale
    "compartment_nippleScale": 5.0,
    # nipple strength
    "compartment_nippleStrength": 10.0,
    # check seeds within radius (mm)
    "compartment_voronSeedRadius": 100.0,

    #####################
    # TDLU variables
    #####################
    # maximum TDLU length
    "tdlu_maxLength": 2.0,
    # minimum TDLU length
    "tdlu_minLength": 1.0,
    # maximum TDLU width
    "tdlu_maxWidth": 1.0,
    # minimum TDLU width
    "tdlu_minWidth": 0.5,

    #####################
    # Perlin noise variables
    #####################
    # maximum fraction of radius deviation
    "perlin_noise_maxDeviation_perlin": 0.1,
    # starting frequency
    "perlin_noise_frequency_perlin": 0.1,
    # octave frequency multiplier
    "perlin_noise_lacunarity_perlin": 2.0,
    # octave signal decay
    "perlin_noise_persistence_perlin": 0.5,
    # number of frequency octaves
    "perlin_noise_numOctaves": 6,
    # x direction noise generation seed
    "perlin_noise_xNoiseGen": 683,
    # y direction noise generation seed
    "perlin_noise_yNoiseGen": 4933,
    # z direction noise generation seed
    "perlin_noise_zNoiseGen": 23,
    # seed noise generation
    "perlin_noise_seedNoiseGen": 3095,
    # shift noise generation seed
    "perlin_noise_shiftNoiseGen": 11,

    #####################
    # Compartment boundary noise
    #####################
    # maximum fraction of distance deviation
    "compartment_boundary_maxDeviation": 0.1,
    # starting frequency
    "compartment_boundary_frequency": 0.15,
    # octave frequency multiplier
    "compartment_boundary_lacunarity": 1.5,
    # octave signal decay
    "compartment_boundary_persistence": 0.5,

    #####################
    # Lobule boundary perturbation noise
    #####################
    # maximum fraction of distance deviation
    "lobule_boundary_maxDeviation": 0.25,
    # starting frequency
    "lobule_boundary_frequency": 0.09,
    # octave frequency multiplier
    "lobule_boundary_lacunarity": 2.0,
    # octave signal decay
    "lobule_boundary_persistence": 0.4,

    #####################
    # Lobule glandular buffer noise
    #####################
    # maximum fraction of distance deviation
    "lobule_glandular_maxDeviation": 0.15,
    # starting frequency
    "lobule_glandular_frequency": 0.05,
    # octave frequency multiplier
    "lobule_glandular_lacunarity": 1.5,
    # octave signal decay
    "lobule_glandular_persistence": 0.5,

    #####################
    # Voronoi segmentation variables
    #####################
    # fat voronoi seed density (mm^-3,)
    "voronoi_fatInFatSeedDensity": 0.001,
    # fat voronoi seed in glandular tissue density (mm^-3,)
    "voronoi_fatInGlandSeedDensity": 0.001,
    # glandular voronoi seed density (mm^-3,)
    "voronoi_glandInGlandSeedDensity": 0.0005,
    # maximum deflection (fraction of pi)
    "voronoi_TDLUDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenTDLU": 0.1,
    # maximum length scale
    "voronoi_maxScaleLenTDLU": 0.2,
    # minimum width scale
    "voronoi_minScaleWidTDLU": 40.0,
    # maximum width scale
    "voronoi_maxScaleWidTDLU": 45.0,
    # minimum strength
    "voronoi_minStrTDLU": 20.0,
    # maximum strength
    "voronoi_maxStrTDLU": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInFatDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInFat": 5.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInFat": 10.0,
    # minimum width scale
    "voronoi_minScaleWidFatInFat": 50.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInFat": 60.0,
    # minimum strength
    "voronoi_minStrFatInFat": 40.0,
    # maximum strength
    "voronoi_maxStrFatInFat": 50.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidFatInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInGland": 40.0,
    # minimum strength
    "voronoi_minStrFatInGland": 20.0,
    # maximum strength
    "voronoi_maxStrFatInGland": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_glandInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenGlandInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenGlandInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidGlandInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidGlandInGland": 40.0,
    # minimum strength
    "voronoi_minStrGlandInGland": 20.0,
    # maximum strength
    "voronoi_maxStrGlandInGland": 22.0,
    # check seeds in radius (mm)
    "voronoi_seedRadius": 40.0,

    #####################
    # fat variables
    #####################
    # min lobule axis length (mm)
    # was 20,
    "fat_minLobuleAxis": 20.0,
    # max lobule axis length (mm)
    # was 35,
    "fat_maxLobuleAxis": 30.0,
    # axial ratio min
    "fat_minAxialRatio_fat": 0.13,
    # axial ratio max
    "fat_maxAxialRatio_fat": 0.75,
    # minimum ligament separation between lobules
    "fat_minLobuleGap": 0.15,
    # maximum of absolute value of Fourier coefficient as fraction of main radius
    "fat_maxCoeffStr": 0.1,
    # minimum of absolute value of Fourier coefficient as fraction of main radius
    "fat_minCoeffStr": 0.05,
    # maximum number of trial lobules
    # was 401,
    "fat_maxLobuleTry": 401,


    #####################
    # ligament variables
    #####################
    "ligament_thickness": 0.1,
    "ligament_targetFrac": 0.85,
    "ligament_maxTry_ligament": 15000,
    "ligament_minAxis": 20.0,
    "ligament_maxAxis": 25.0,
    "ligament_minAxialRatio_ligament": 0.2,
    "ligament_maxAxialRatio_ligament": 0.3,
    "ligament_maxPerturb": 0.05,
    "ligament_maxDeflect": 0.12,
    "ligament_scale": 0.007,
    "ligament_lacunarity": 1.5,
    "ligament_persistence": 0.3,
    "ligament_numOctaves": 6,

    #####################
    # duct tree variables
    #####################
    # target number of branches (uint)
    "duct_tree_maxBranch": 400,
    # maximum generation (uint)
    "duct_tree_maxGen": 7,
    # initial radius of tree (mm)
    "duct_tree_initRad": 0.5,
    # base Length of root duct at nipple (mm)
    "duct_tree_baseLength": 7.6,
    # number of voxels for tree density tracking (uint)
    "duct_tree_nFillX": 50,
    "duct_tree_nFillY": 50,
    "duct_tree_nFillZ": 50,

    #####################
    # duct branch variables
    #####################
    # minimum branch radius to have children (mm)
    "duct_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "duct_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "duct_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "duct_branch_lenShrink": 0.5,
    # maximum jitter in branch length (fraction)
    "duct_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "duct_branch_rotateJitter": 0.1,

    #####################
    # duct segment variables
    #####################
    # radius distribution shape parameters
    "duct_segment_radiusBetaA": 6.0,
    "duct_segment_radiusBetaB": 10.0,
    # fraction of branch length per segment
    "duct_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "duct_segment_maxCurvRad": 10.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "duct_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "duct_segment_minEndRad": 0.95,
    "duct_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "duct_segment_angleWt": 1.0,
    # cost function density weighting
    "duct_segment_densityWt": 20.0,
    # number of trial segments to generate (uint)
    "duct_segment_numTry": 50,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "duct_segment_maxTry": 100,
    # total number of segment tries before completely giving up
    "duct_segment_absMaxTry": 10000,
    # step size for checking segment is valid (mm)
    "duct_segment_roiStep": 0.1,


    #####################
    # vessel tree variables
    #####################
    # target number of branches (uint)
    "vessel_tree_maxBranch": 750,
    # maximum generation (uint)
    "vessel_tree_maxGen": 6,
    # initial radius of tree (mm)
    "vessel_tree_initRad": 0.75,
    # base length of root vessel (mm)
    "vessel_tree_baseLength": 15.0,
    # number of voxels for tree density tracking (uint)
    "vessel_tree_nFillX": 30,
    "vessel_tree_nFillY": 69,
    "vessel_tree_nFillZ": 69,

    #####################
    # vessel branch variables
    #####################
    # minimum branch radius to have children (mm)
    "vessel_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "vessel_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "vessel_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "vessel_branch_lenShrink": 0.8,
    # maximum jitter in branch length (fraction)
    "vessel_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "vessel_branch_rotateJitter": 0.1,

    #####################
    # vessel segment variables
    #####################
    # radius distribution shape parameters
    "vessel_segment_radiusBetaA": 6.0,
    "vessel_segment_radiusBetaB": 10.0,
    # fraction of branch length to segment
    "vessel_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "vessel_segment_maxCurvRad": 200.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "vessel_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "vessel_segment_minEndRad": 0.95,
    "vessel_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "vessel_segment_angleWt": 100.0,
    # cost function density weighting
    "vessel_segment_densityWt": 1.0,
    # cost function direction weighting
    "vessel_segment_dirWt": 100.0,
    # number of trial segments to generate (uint)
    "vessel_segment_numTry": 100,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "vessel_segment_maxTry": 300,
    # total number of segment tries before completely giving up
    "vessel_segment_absMaxTry": 100000,
    # step size for checking segment is valid (mm)
    "vessel_segment_roiStep": 0.1,
}

VICTRE_FATTY = {
    "compressionThickness": 60,  # mm
    # phantom voxel size (mm)
    "imgRes": 0.05,
    # thickness of breast skin (mm)
    "skinThick": 0.75,
    # nipple length (mm)
    "nippleLen": 4.0,
    # nipple radius (mm)
    "nippleRad": 4.0,
    # nipple radius (mm)
    "areolaRad": 8.0,
    # left breast - select left or right breast (boolean)
    "leftBreast": "false",
    # desired fat fraction
    "targetFatFrac": 0.95,

    #####################
    # breast surface shape
    #####################
    # u resolution of base shape
    "surface_ures": 0.005,
    # v resolution of base shape
    "surface_vres": 0.005,
    # minimum point separation (mm)
    "surface_pointSep": 0.005,
    # back ring thickness (mm)
    "surface_ringWidth": 10.0,
    # back ring step size (mm)
    "surface_ringSep": 0.5,
    # angle to preserve while smoothing (degrees)
    "surface_featureAngle": 20.0,
    # fraction of triangles to decimate
    "surface_targetReduction": 0.05,
    # bottom scale was 1.3,
    "surface_a1b": 1.4,
    # top scale was 1.3,
    "surface_a1t": 1.4,
    # left scale
    "surface_a2l": 1.4,
    # right scale
    "surface_a2r": 1.4,
    # outward scale was 1.75,
    "surface_a3": 2.75,
    # u quadric exponent
    "surface_eps1": 1.1,
    # v quadric exponent
    "surface_eps2": 1.0,
    # do ptosis deformation (boolean)
    "surface_doPtosis": "true",
    "surface_ptosisB0": 0.2,
    "surface_ptosisB1": 0.05,
    # do turn deformation (boolean)
    "surface_doTurn": "false",
    "surface_turnC0": -0.498,
    "surface_turnC1": 0.213,
    # do top shape deformation (boolean)
    "surface_doTopShape": "true",
    "surface_topShapeS0": 0.0,
    "surface_topShapeS1": 0.0,
    "surface_topShapeT0": -8.0,
    "surface_topShapeT1": -2.0,
    # do flatten size deformation (boolean)
    "surface_doFlattenSide": "true",
    "surface_flattenSideG0": 1.5,
    "surface_flattenSideG1": -0.5,
    # do turn top deformation (boolean)
    "surface_doTurnTop": "true",
    "surface_turnTopH0": 0.166,
    "surface_turnTopH1": -0.372,

    #####################
    # breast compartment
    #####################
    # number of breast compartments
    "compartment_num": 10,
    # distance along nipple line of compartment seed base (mm)
    "compartment_seedBaseDist": 16,
    # fraction of phantom in nipple direction forced to be fat
    "compartment_backFatBufferFrac": 0.01,
    # number of backplane seed points
    "compartment_numBackSeeds": 250,
    # maximum seed jitter (fraction of subtended angle)
    "compartment_angularJitter": 0.125,
    # maximum seed jitter in nipple direction (mm)
    "compartment_zJitter": 5.0,
    # maximum radial distance from base seed as a fraction of distance to breast surface
    "compartment_maxFracRadialDist": 0.5,
    # minimum radial distance from base seed as a fraction of distance to breast surface
    "compartment_minFracRadialDist": 0.25,
    # minimum scale in nipple direction
    "compartment_minScaleNippleDir": 0.01,
    # maximum scale in nipple direction
    "compartment_maxScaleNippleDir": 0.01,
    # minimum scale in non-nipple direction
    "compartment_minScale": 30.0,
    # maximum scale in non-nipple direction
    "compartment_maxScale": 40.0,
    # minimum gland strength
    "compartment_minGlandStrength": 30.0,
    # maximum gland strength
    "compartment_maxGlandStrength": 30.0,
    # maximum compartment deflection angle from pointing towards nipple (fraction of pi)
    "compartment_maxDeflect": 0.01,
    # minimum scale skin seeds in nipple direction
    "compartment_minSkinScaleNippleDir": 5.0,
    # maximum scale skin seeds in nipple direction
    "compartment_maxSkinScaleNippleDir": 5.0,
    # minimum scale skin in non-nipple direction
    "compartment_minSkinScale": 200.0,
    # maximum scale skin in non-nipple direction
    "compartment_maxSkinScale": 200.0,
    # skin strength
    "compartment_skinStrength": 2.0,
    # back scale
    "compartment_backScale": 60.0,
    # back strength
    "compartment_backStrength": 4.0,
    # nipple scale
    "compartment_nippleScale": 5.0,
    # nipple strength
    "compartment_nippleStrength": 10.0,
    # check seeds within radius (mm)
    "compartment_voronSeedRadius": 100.0,

    #####################
    # TDLU variables
    #####################
    # maximum TDLU length
    "tdlu_maxLength": 2.0,
    # minimum TDLU length
    "tdlu_minLength": 1.0,
    # maximum TDLU width
    "tdlu_maxWidth": 1.0,
    # minimum TDLU width
    "tdlu_minWidth": 0.5,

    #####################
    # Perlin noise variables
    #####################
    # maximum fraction of radius deviation
    "perlin_noise_maxDeviation": 0.1,
    # starting frequency
    "perlin_noise_frequency": 0.1,
    # octave frequency multiplier
    "perlin_noise_lacunarity": 2.0,
    # octave signal decay
    "perlin_noise_persistence": 0.5,
    # number of frequency octaves
    "perlin_noise_numOctaves": 6,
    # x direction noise generation seed
    "perlin_noise_xNoiseGen": 683,
    # y direction noise generation seed
    "perlin_noise_yNoiseGen": 4933,
    # z direction noise generation seed
    "perlin_noise_zNoiseGen": 23,
    # seed noise generation
    "perlin_noise_seedNoiseGen": 3095,
    # shift noise generation seed
    "perlin_noise_shiftNoiseGen": 11,

    #####################
    # Compartment boundary noise
    #####################
    # maximum fraction of distance deviation
    "compartment_boundary_maxDeviation": 0.1,
    # starting frequency
    "compartment_boundary_frequency": 0.15,
    # octave frequency multiplier
    "compartment_boundary_lacunarity": 1.5,
    # octave signal decay
    "compartment_boundary_persistence": 0.5,

    #####################
    # Lobule boundary perturbation noise
    #####################
    # maximum fraction of distance deviation
    "lobule_boundary_maxDeviation": 0.25,
    # starting frequency
    "lobule_boundary_frequency": 0.09,
    # octave frequency multiplier
    "lobule_boundary_lacunarity": 2.0,
    # octave signal decay
    "lobule_boundary_persistence": 0.4,

    #####################
    # Lobule glandular buffer noise
    #####################
    # maximum fraction of distance deviation
    "lobule_glandular_maxDeviation": 0.15,
    # starting frequency
    "lobule_glandular_frequency": 0.05,
    # octave frequency multiplier
    "lobule_glandular_lacunarity": 1.5,
    # octave signal decay
    "lobule_glandular_persistence": 0.5,

    #####################
    # Voronoi segmentation variables
    #####################
    # fat voronoi seed density (mm^-3,)
    "voronoi_fatInFatSeedDensity": 0.001,
    # fat voronoi seed in glandular tissue density (mm^-3,)
    "voronoi_fatInGlandSeedDensity": 0.001,
    # glandular voronoi seed density (mm^-3,)
    "voronoi_glandInGlandSeedDensity": 0.0005,
    # maximum deflection (fraction of pi)
    "voronoi_TDLUDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenTDLU": 0.1,
    # maximum length scale
    "voronoi_maxScaleLenTDLU": 0.2,
    # minimum width scale
    "voronoi_minScaleWidTDLU": 40.0,
    # maximum width scale
    "voronoi_maxScaleWidTDLU": 45.0,
    # minimum strength
    "voronoi_minStrTDLU": 20.0,
    # maximum strength
    "voronoi_maxStrTDLU": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInFatDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInFat": 5.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInFat": 10.0,
    # minimum width scale
    "voronoi_minScaleWidFatInFat": 50.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInFat": 60.0,
    # minimum strength
    "voronoi_minStrFatInFat": 40.0,
    # maximum strength
    "voronoi_maxStrFatInFat": 50.0,
    # maximum deflection (fraction of pi)
    "voronoi_fatInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenFatInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenFatInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidFatInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidFatInGland": 40.0,
    # minimum strength
    "voronoi_minStrFatInGland": 20.0,
    # maximum strength
    "voronoi_maxStrFatInGland": 22.0,
    # maximum deflection (fraction of pi)
    "voronoi_glandInGlandDeflectMax": 0.15,
    # minimum length scale
    "voronoi_minScaleLenGlandInGland": 1.0,
    # maximum length scale
    "voronoi_maxScaleLenGlandInGland": 2.0,
    # minimum width scale
    "voronoi_minScaleWidGlandInGland": 30.0,
    # maximum width scale
    "voronoi_maxScaleWidGlandInGland": 40.0,
    # minimum strength
    "voronoi_minStrGlandInGland": 20.0,
    # maximum strength
    "voronoi_maxStrGlandInGland": 22.0,
    # check seeds in radius (mm)
    "voronoi_seedRadius": 40.0,

    #####################
    # fat variables
    #####################
    # min lobule axis length (mm)
    # was 20,
    "fat_minLobuleAxis": 20.0,
    # max lobule axis length (mm)
    # was 35,
    "fat_maxLobuleAxis": 30.0,
    # axial ratio min
    "fat_minAxialRatio": 0.13,
    # axial ratio max
    "fat_maxAxialRatio": 0.75,
    # minimum ligament separation between lobules
    "fat_minLobuleGap": 0.15,
    # maximum of absolute value of Fourier coefficient as fraction of main radius
    "fat_maxCoeffStr": 0.1,
    # minimum of absolute value of Fourier coefficient as fraction of main radius
    "fat_minCoeffStr": 0.05,
    # maximum number of trial lobules
    # was 401,
    "fat_maxLobuleTry": 401,


    #####################
    # ligament variables
    #####################
    "ligament_thickness": 0.1,
    "ligament_targetFrac": 0.85,
    "ligament_maxTry": 15000,
    "ligament_minAxis": 20.0,
    "ligament_maxAxis": 25.0,
    "ligament_minAxialRatio": 0.2,
    "ligament_maxAxialRatio": 0.3,
    "ligament_maxPerturb": 0.05,
    "ligament_maxDeflect": 0.12,
    "ligament_scale": 0.007,
    "ligament_lacunarity": 1.5,
    "ligament_persistence": 0.3,
    "ligament_numOctaves": 6,

    #####################
    # duct tree variables
    #####################
    # target number of branches (uint)
    "duct_tree_maxBranch": 400,
    # maximum generation (uint)
    "duct_tree_maxGen": 7,
    # initial radius of tree (mm)
    "duct_tree_initRad": 0.5,
    # base Length of root duct at nipple (mm)
    "duct_tree_baseLength": 19.0,
    # number of voxels for tree density tracking (uint)
    "duct_tree_nFillX": 50,
    "duct_tree_nFillY": 50,
    "duct_tree_nFillZ": 50,

    #####################
    # duct branch variables
    #####################
    # minimum branch radius to have children (mm)
    "duct_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "duct_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "duct_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "duct_branch_lenShrink": 0.5,
    # maximum jitter in branch length (fraction)
    "duct_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "duct_branch_rotateJitter": 0.1,

    #####################
    # duct segment variables
    #####################
    # radius distribution shape parameters
    "duct_segment_radiusBetaA": 6.0,
    "duct_segment_radiusBetaB": 10.0,
    # fraction of branch length per segment
    "duct_segment_segFrac": 0.15,
    # maximum radius of curvature (mm)
    "duct_segment_maxCurvRad": 10.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "duct_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "duct_segment_minEndRad": 0.95,
    "duct_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "duct_segment_angleWt": 1.0,
    # cost function density weighting
    "duct_segment_densityWt": 20.0,
    # number of trial segments to generate (uint)
    "duct_segment_numTry": 50,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "duct_segment_maxTry": 100,
    # total number of segment tries before completely giving up
    "duct_segment_absMaxTry": 10000,
    # step size for checking segment is valid (mm)
    "duct_segment_roiStep": 0.1,


    #####################
    # vessel tree variables
    #####################
    # target number of branches (uint)
    "vessel_tree_maxBranch": 750,
    # maximum generation (uint)
    "vessel_tree_maxGen": 6,
    # initial radius of tree (mm)
    "vessel_tree_initRad": 0.75,
    # base length of root vessel (mm)
    "vessel_tree_baseLength": 15.0,
    # number of voxels for tree density tracking (uint)
    "vessel_tree_nFillX": 30,
    "vessel_tree_nFillY": 69,
    "vessel_tree_nFillZ": 69,

    #####################
    # vessel branch variables
    #####################
    # minimum branch radius to have children (mm)
    "vessel_branch_childMinRad": 0.1,
    # minimum starting radius as a fraction of parent end radius
    "vessel_branch_minRadFrac": 0.65,
    # maximum starting radius as a fraction of parent end radius
    "vessel_branch_maxRadFrac": 0.99,
    # length reduction as fraction of parent length
    "vessel_branch_lenShrink": 0.8,
    # maximum jitter in branch length (fraction)
    "vessel_branch_lenRange": 0.1,
    # aximuthal angle noise (radians)
    "vessel_branch_rotateJitter": 0.1,

    #####################
    # vessel segment variables
    #####################
    # radius distribution shape parameters
    "vessel_segment_radiusBetaA": 6.0,
    "vessel_segment_radiusBetaB": 10.0,
    # fraction of branch length to segment
    "vessel_segment_segFrac": 0.25,
    # maximum radius of curvature (mm)
    "vessel_segment_maxCurvRad": 200.0,
    # maximum length of segment based on
    # curvature (fraction of pi radians)
    "vessel_segment_maxCurvFrac": 0.5,
    # min and max end radius as fraction of start radius
    "vessel_segment_minEndRad": 0.95,
    "vessel_segment_maxEndRad": 1.0,
    # cost function preferential angle weighting
    "vessel_segment_angleWt": 100.0,
    # cost function density weighting
    "vessel_segment_densityWt": 1.0,
    # cost function direction weighting
    "vessel_segment_dirWt": 100.0,
    # number of trial segments to generate (uint)
    "vessel_segment_numTry": 100,
    # maximum number of segments to generate before
    # giving up and reducing length (uint)
    "vessel_segment_maxTry": 300,
    # total number of segment tries before completely giving up
    "vessel_segment_absMaxTry": 100000,
    # step size for checking segment is valid (mm)
    "vessel_segment_roiStep": 0.1,
}

DENSITY_RANGES = {
    "breastHeight":[71.7,89.7,108.4,109.9],
    "targetFatFrac": [0.4, 0.66, 0.85, 0.95],
    "surface_a1b": [1.0, 1.2, 1.4, 1.4],
    "surface_a1t": [1.0, 1.2, 1.4, 1.4],
    "surface_a2l": [1.0, 1.0, 1.2, 1.4],
    "surface_a2r": [1.0, 1.0, 1.2, 1.4],
    "surface_a3": [1.0, 1.55, 2.0, 2.75],
    "surface_eps1": [1.2, 1.2, 1.2, 1.1],
    "surface_topShapeT0": [-12.0, -12.0, -8.0, -8.0],
    "surface_topShapeT1": [-5.0, -5.0, -2.0, -2.0],
    "compartment_backFatBufferFrac": [0.008, 0.01, 0.01, 0.01],
    "compartment_numBackSeeds": [150, 250, 250, 250],
    "compartment_minSkinScaleNippleDir": [10.0, 10.0, 5.0, 5.0],
    "compartment_maxSkinScaleNippleDir": [20.0, 20.0, 5.0, 5.0],
    "compartment_maxSkinScale":[400,200,200,200],
    "compartment_skinStrength":[0.5,1.0,2.0,2.0],
    "compartment_backStrength":[1.0,2.0,4.0,4.0],
    "duct_tree_baseLength":[7.6,7.6,7.6,19],
    "duct_segment_segFrac":[0.25,0.25,0.25,0.15],
    # mm (this is the objective thickness)
    "compressionThickness": [35, 45, 55, 60],
    "number_histories": [7.8e9, 1.02e10, 2.04e10, 2.22e10],
}

MHD_FILE = """ObjectType = $ObjectType
NDims = $NDims
BinaryData = $BinaryData
BinaryDataByteOrderMSB = $BinaryDataByteOrderMSB
CompressedData = $CompressedData
TransformMatrix = $TransformMatrix
Offset = $Offset
CenterOfRotation = $CenterOfRotation
ElementSpacing = $ElementSpacing
DimSize = $DimSize
AnatomicalOrientation = $AnatomicalOrientation
ElementType = $ElementType
ObjectType = $ObjectType
ElementDataFile = $ElementDataFile"""

# mass configuration file
VICTRE_DEFAULT_SPICULATED_MASS = {
    # [base]
    # voxel size (mm)
    "imgRes": 0.05,
    # complexity scaling for VICTRE 0:simple, 1:full complexity
    "complexity": 1.0,
    # random number seed (unsigned int) from /dev/urandom if unspecified
    "seed": 0,

    # [mass]
    # maximum spherical harmonic order
    "lMax": 4,
    # mean mass radius (mm)
    "alpha": 2.5,
    # mean mass surface irregularity variance (mm^2)
    "meanSigma2": 0.31,
    # mass surface irregularity standard deviation (mm^2)
    "stdSigma2": 0.04,
    # covariance power law index
    "powerLaw": 4.0,
    # mean number of low freq. modifications
    "meanLF": 611.2,
    # std. deviation of number of low freq. modifications
    "stdLF": 70.6,
    # mean LF shape distribution 0:spike, 1:bump
    "meanShape": 0.36,
    # std. deviation of LF shape 0:spike, 1:bump
    "stdShape": 0.48,
    # mean relative LF radius
    "meanLFRad": 0.229,
    # std. deviation of LF radius
    "stdLFRad": 0.073,
    # mean relative LF length
    "meanLFLen": 0.113,
    # std. deviation of LF length
    "stdLFLen": 0.021,
    # mean fuzzy alpha
    "meanFuzzAlpha": 0.015,
    # std. deviation of fuzzy alpha
    "stdFuzzAlpha": 0.0,

    # [spicule]
    # mean number of initial segments
    "meanInitial": 1358,
    # std. deviation number of initial segments
    "stdInitial": 365,
    # mean max number of neighbor segments
    "meanNeigh": 8.98,
    # std. deviation max number of neighbor segments
    "stdNeigh": 1.89,
    # mean initial relative radius
    "meanInitRad": 0.024,
    # std. deviation initial relative radius
    "stdInitRad": 0.0053,
    # mean radius decrease
    "meanRadDec": 0.89,
    # std. deviation radius decrease
    "stdRadDec": 0.31,
    # mean initial relative length
    "meanInitLen": 0.173,
    # std. deviation initial relative length
    "stdInitLen": 0.018,
    # mean length decrease
    "meanLenDec": 0.91,
    # std. deviation length decrease
    "stdLenDec": 0.30,
    # mean continue prob.
    "meanContProb": 0.717,
    # std. deviation continue prob.
    "stdContProb": 0.057,
    # mean branch angle.
    "meanBranchAng": 6.55,
    # std. deviation branch angle.
    "stdBranchAng": 0.62
}

DICOM_SCALING={"Siemens":
    {
        "meanAdditiveNoise":5200,
        "conversionFactorDM":0.000239,
        "conversionFactorDBT":0.001469,
        "offset":50,
        "toUInt16":8
    }
}