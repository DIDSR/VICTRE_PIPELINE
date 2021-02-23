"
# 
# >>>> INPUT FILE FOR MC-GPU v1.5 VICTRE-DBT >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#                
#  This input file simulates a mammogram and 25 projections of a DBT scan (+-25deg). 
#  Main acquistion parameters: 
#     - Source-to-detector distance 65 cm. 
#     - Pixel size 85 micron (= 25.5 cm / 3000 pixels)
#     - Antiscatter grid used only in the mammogram; motion blur used only in the DBT scan.
#     - Breast phantom must be generated using C. Graff's software (hardcoded conversion from binary voxel value to material and density)
#     - Number of histories matches number of x rays in a DBT projection, to reproduce the quantum noise and dose.
#        -- Mammogram simulated with 2/3 the histories in the 25 projections combined (1.02e10*25*2/3 histories)
#        -- It is ok to reduce the number of histories for testing the code!
#        -- Number of histories computed to match the air kerma measured with the real system at the center of a PMMA phantom of equivalent thickness.
#        -- Number of histotries must be re-calculated if the energy spectrum or beam aperture (field size) are changed.
#
#                      [Andreu Badal, 2019-08-23]
#

#[SECTION SIMULATION CONFIG v.2009-05-12]
$number_histories                         # TOTAL NUMBER OF HISTORIES, OR SIMULATION TIME IN SECONDS IF VALUE < 100000
$random_seed                        # RANDOM SEED (ranecu PRNG)
$selected_gpu                              # GPU NUMBER TO USE WHEN MPI IS NOT USED, OR TO BE AVOIDED IN MPI RUNS
$gpu_threads                             # GPU THREADS PER CUDA BLOCK (multiple of 32)
$histories_per_thread                            # SIMULATED HISTORIES PER GPU THREAD
 
#[SECTION SOURCE v.2016-12-02]
$spectrum_file # X-RAY ENERGY SPECTRUM FILE
$source_position           # SOURCE POSITION: X (chest-to-nipple), Y (right-to-left), Z (caudal-to-cranial) [cm]
$source_direction             # SOURCE DIRECTION COSINES: U V W
$fam_beam_aperture    # ==> 2/3 original angle of 11.203       # TOTAL AZIMUTHAL (WIDTH, X) AND POLAR (HEIGHT, Z) APERTURES OF THE FAN BEAM [degrees] (input negative to automatically cover the whole detector)
$euler_angles             # EULER ANGLES (RzRyRz) TO ROTATE RECTANGULAR BEAM FROM DEFAULT POSITION AT Y=0, NORMAL=(0,-1,0)
$focal_spot                         # SOURCE GAUSSIAN FOCAL SPOT FWHM [cm]
$angular_blur                           # 0.18 for DBT, 0 for FFDM [Mackenzie2017]  # ANGULAR BLUR DUE TO MOVEMENT ([exposure_time]*[angular_speed]) [degrees]
$collimate_beam                             # COLLIMATE BEAM TOWARDS POSITIVE AZIMUTHAL (X) ANGLES ONLY? (ie, cone-beam center aligned with chest wall in mammography) [YES/NO]
 
#[SECTION IMAGE DETECTOR v.2017-06-20]
$output_file   # OUTPUT IMAGE FILE NAME
$image_pixels                  # NUMBER OF PIXELS IN THE IMAGE: Nx Nz
$image_size                 # IMAGE SIZE (width, height): Dx Dz [cm]
$distance_source                           # SOURCE-TO-DETECTOR DISTANCE (detector set in front of the source, perpendicular to the initial direction)
$image_offset                     # IMAGE OFFSET ON DETECTOR PLANE IN WIDTH AND HEIGHT DIRECTIONS (BY DEFAULT BEAM CENTERED AT IMAGE CENTER) [cm]
$detector_thickness                         # DETECTOR THICKNESS [cm]
$mean_free_path  # ==> MFP(Se,19.0keV)   # DETECTOR MATERIAL MEAN FREE PATH AT AVERAGE ENERGY [cm]
$k_edge_energy  # DETECTOR K-EDGE ENERGY [eV], K-FLUORESCENCE ENERGY [eV], K-FLUORESCENCE YIELD, MFP AT FLUORESCENCE ENERGY [cm]
$detector_gain                   # EFECTIVE DETECTOR GAIN, W_+- [eV/ehp], AND SWANK FACTOR (input 0 to report ideal energy fluence)
$additive_noise                         # ADDITIVE ELECTRONIC NOISE LEVEL (electrons/pixel)
$cover_thickness          # ==> MFP(polystyrene,19keV)       # PROTECTIVE COVER THICKNESS (detector+grid) [cm], MEAN FREE PATH AT AVERAGE ENERGY [cm]
$antiscatter_grid_ratio            # ANTISCATTER GRID RATIO, FREQUENCY, STRIP THICKNESS [X:1, lp/cm, cm] (enter 0 to disable the grid)
$antiscatter_strips   # ==> MFP(lead&polystyrene,19keV)  # ANTISCATTER STRIPS AND INTERSPACE MEAN FREE PATHS AT AVERAGE ENERGY [cm]
$antiscatter_grid_lines                              # ORIENTATION 1D FOCUSED ANTISCATTER GRID LINES: 0==STRIPS PERPENDICULAR LATERAL DIRECTION (mammo style); 1==STRIPS PARALLEL LATERAL DIRECTION (DBT style)

#[SECTION TOMOGRAPHIC TRAJECTORY v.2016-12-02]
$number_projections      # ==> 1 for mammo only; ==> 25 for mammo + DBT    # NUMBER OF PROJECTIONS (1 disables the tomographic mode)
$rotation_axis_distance                            # SOURCE-TO-ROTATION AXIS DISTANCE
$projections_angle           # ANGLE BETWEEN PROJECTIONS (360/num_projections for full CT) [degrees]
$angular_rotation_first                           # ANGULAR ROTATION TO FIRST PROJECTION (USEFUL FOR DBT, INPUT SOURCE DIRECTION CONSIDERED AS 0 DEGREES) [degrees]
$rotation_axis                  # AXIS OF ROTATION (Vx,Vy,Vz)
$axis_translation                            # TRANSLATION ALONG ROTATION AXIS BETWEEN PROJECTIONS (HELICAL SCAN) [cm]
$detector_fixed                             # KEEP DETECTOR FIXED AT 0 DEGREES FOR DBT? [YES/NO]
$simulate_both                             # SIMULATE BOTH 0 deg PROJECTION AND TOMOGRAPHIC SCAN (WITHOUT GRID) WITH 2/3 TOTAL NUM HIST IN 1st PROJ (eg, DBT+mammo)? [YES/NO]

#[SECTION DOSE DEPOSITION v.2012-12-12]
$tally_material_dose                             # TALLY MATERIAL DOSE? [YES/NO] (electrons not transported, x-ray energy locally deposited at interaction)
$tally_voxel_dose                              # TALLY 3D VOXEL DOSE? [YES/NO] (dose measured separately for each voxel)
$output_dose_filename                 # OUTPUT VOXEL DOSE FILE NAME
$roi_voxel_dose_x                        # VOXEL DOSE ROI: X-index min max (first voxel has index 1)
$roi_voxel_dose_y                        # VOXEL DOSE ROI: Y-index min max
$roi_voxel_dose_z                        # VOXEL DOSE ROI: Z-index min max
 
#[SECTION VOXELIZED GEOMETRY FILE v.2017-07-26]
$phantom_file    # VOXEL GEOMETRY FILE (penEasy 2008 format; .gz accepted)
$voxel_geometry_offset              # OFFSET OF THE VOXEL GEOMETRY (DEFAULT ORIGIN AT LOWER BACK CORNER) [cm]
$number_voxels              # NUMBER OF VOXELS: INPUT A 0 TO READ ASCII FORMAT WITH HEADER SECTION, RAW VOXELS WILL BE READ OTHERWISE
$voxel_size           # VOXEL SIZES [cm]
$low_resolution_voxel_size                          # SIZE OF LOW RESOLUTION VOXELS THAT WILL BE DESCRIBED BY A BINARY TREE, GIVEN AS POWERS OF TWO (eg, 2 2 3 = 2^2x2^2x2^3 = 128 input voxels per low res voxel; 0 0 0 disables tree)
 
#[SECTION MATERIAL FILE LIST v.2009-11-30]
