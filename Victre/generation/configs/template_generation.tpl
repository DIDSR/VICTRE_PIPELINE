##########################
# breast phantom configuration file
##########################


#####################
# basic variables
#####################
[base]
# output directory
outputDir=$outputDir
# phantom voxel size (mm)
imgRes=$imgRes
# thickness of breast skin (mm)
skinThick=$skinThick
# nipple length (mm)
nippleLen=$nippleLen
# nipple radius (mm)
nippleRad=$nippleRad
# nipple radius (mm)
areolaRad=$areolaRad
# left breast - select left or right breast (boolean)
leftBreast=$leftBreast
# desired fat fraction
targetFatFrac=$targetFatFrac
# random number seed (unsigned int)
# chosen randomly if not set
seed=$seed

#####################
# breast surface shape
#####################
[shape]
# u resolution of base shape
ures=$surface_ures
# v resolution of base shape
vres=$surface_vres
# minimum point separation (mm)
pointSep=$surface_pointSep
# back ring thickness (mm)
ringWidth=$surface_ringWidth
# back ring step size (mm)
ringSep=$surface_ringSep
# angle to preserve while smoothing (degrees)
featureAngle=$surface_featureAngle
# fraction of triangles to decimate
targetReduction=$surface_targetReduction
# bottom scale was 1.3
a1b=$surface_a1b
# top scale was 1.3
a1t=$surface_a1t
# left scale
a2l=$surface_a2l
# right scale
a2r=$surface_a2r
# outward scale was 1.75
a3=$surface_a3
# u quadric exponent
eps1=$surface_eps1
# v quadric exponent
eps2=$surface_eps2
# do ptosis deformation (boolean)
doPtosis=$surface_doPtosis
ptosisB0=$surface_ptosisB0
ptosisB1=$surface_ptosisB1
# do turn deformation (boolean)
doTurn=$surface_doTurn
turnC0=$surface_turnC0
turnC1=$surface_turnC1
# do top shape deformation (boolean)
doTopShape=$surface_doTopShape
topShapeS0=$surface_topShapeS0
topShapeS1=$surface_topShapeS1
topShapeT0=$surface_topShapeT0
topShapeT1=$surface_topShapeT1
# do flatten size deformation (boolean)
doFlattenSide=$surface_doFlattenSide
flattenSideG0=$surface_flattenSideG0
flattenSideG1=$surface_flattenSideG1
# do turn top deformation (boolean)
doTurnTop=$surface_doTurnTop
turnTopH0=$surface_turnTopH0
turnTopH1=$surface_turnTopH1

#####################
# breast compartment
#####################
[compartments]
# number of breast compartments
num=$compartment_num
# distance along nipple line of compartment seed base (mm)
seedBaseDist=$compartment_seedBaseDist
# fraction of phantom in nipple direction forced to be fat
backFatBufferFrac=$compartment_backFatBufferFrac
# number of backplane seed points
numBackSeeds=$compartment_numBackSeeds
# maximum seed jitter (fraction of subtended angle)
angularJitter=$compartment_angularJitter
# maximum seed jitter in nipple direction (mm)
zJitter=$compartment_zJitter
# maximum radial distance from base seed as a fraction of distance to breast surface
maxFracRadialDist=$compartment_maxFracRadialDist
# minimum radial distance from base seed as a fraction of distance to breast surface
minFracRadialDist=$compartment_minFracRadialDist
# minimum scale in nipple direction
minScaleNippleDir=$compartment_minScaleNippleDir
# maximum scale in nipple direction
maxScaleNippleDir=$compartment_maxScaleNippleDir
# minimum scale in non-nipple direction
minScale=$compartment_minScale
# maximum scale in non-nipple direction
maxScale=$compartment_maxScale
# minimum gland strength
minGlandStrength=$compartment_minGlandStrength
# maximum gland strength
maxGlandStrength=$compartment_maxGlandStrength
# maximum compartment deflection angle from pointing towards nipple (fraction of pi)
maxDeflect=$compartment_maxDeflect
# minimum scale skin seeds in nipple direction
minSkinScaleNippleDir=$compartment_minSkinScaleNippleDir
# maximum scale skin seeds in nipple direction
maxSkinScaleNippleDir=$compartment_maxSkinScaleNippleDir
# minimum scale skin in non-nipple direction
minSkinScale=$compartment_minSkinScale
# maximum scale skin in non-nipple direction
maxSkinScale=$compartment_maxSkinScale
# skin strength
skinStrength=$compartment_skinStrength
# back scale
backScale=$compartment_backScale
# back strength
backStrength=$compartment_backStrength
# nipple scale
nippleScale=$compartment_nippleScale
# nipple strength
nippleStrength=$compartment_nippleStrength
# check seeds within radius (mm)
voronSeedRadius=$compartment_voronSeedRadius

#####################
# TDLU variables
#####################
[TDLU]
# maximum TDLU length
maxLength=$tdlu_maxLength
# minimum TDLU length
minLength=$tdlu_minLength
# maximum TDLU width
maxWidth=$tdlu_maxWidth
# minimum TDLU width
minWidth=$tdlu_minWidth

#####################
# Perlin noise variables
#####################
[perlin]
# maximum fraction of radius deviation 
maxDeviation=$perlin_noise_maxDeviation
# starting frequency
frequency=$perlin_noise_frequency
# octave frequency multiplier
lacunarity=$perlin_noise_lacunarity
# octave signal decay
persistence=$perlin_noise_persistence
# number of frequency octaves
numOctaves=$perlin_noise_numOctaves
# x direction noise generation seed
xNoiseGen=$perlin_noise_xNoiseGen
# y direction noise generation seed
yNoiseGen=$perlin_noise_yNoiseGen
# z direction noise generation seed
zNoiseGen=$perlin_noise_zNoiseGen
# seed noise generation
seedNoiseGen=$perlin_noise_seedNoiseGen
# shift noise generation seed
shiftNoiseGen=$perlin_noise_shiftNoiseGen

#####################
# Compartment boundary noise
#####################
[boundary]
# maximum fraction of distance deviation 
maxDeviation=$compartment_boundary_maxDeviation
# starting frequency
frequency=$compartment_boundary_frequency
# octave frequency multiplier
lacunarity=$compartment_boundary_lacunarity
# octave signal decay
persistence=$compartment_boundary_persistence

#####################
# Lobule boundary perturbation noise
#####################
[perturb]
# maximum fraction of distance deviation 
maxDeviation=$lobule_boundary_maxDeviation
# starting frequency
frequency=$lobule_boundary_frequency
# octave frequency multiplier
lacunarity=$lobule_boundary_lacunarity
# octave signal decay
persistence=$lobule_boundary_persistence

#####################
# Lobule glandular buffer noise
#####################
[buffer]
# maximum fraction of distance deviation 
maxDeviation=$lobule_glandular_maxDeviation
# starting frequency
frequency=$lobule_glandular_frequency
# octave frequency multiplier
lacunarity=$lobule_glandular_lacunarity
# octave signal decay
persistence=$lobule_glandular_persistence

#####################
# Voronoi segmentation variables
#####################
[voronoi]
# fat voronoi seed density (mm^-3)
fatInFatSeedDensity=$voronoi_fatInFatSeedDensity
# fat voronoi seed in glandular tissue density (mm^-3)
fatInGlandSeedDensity=$voronoi_fatInGlandSeedDensity
# glandular voronoi seed density (mm^-3)
glandInGlandSeedDensity=$voronoi_glandInGlandSeedDensity
# maximum deflection (fraction of pi)
TDLUDeflectMax=$voronoi_TDLUDeflectMax
# minimum length scale
minScaleLenTDLU=$voronoi_minScaleLenTDLU
# maximum length scale
maxScaleLenTDLU=$voronoi_maxScaleLenTDLU
# minimum width scale
minScaleWidTDLU=$voronoi_minScaleWidTDLU
# maximum width scale
maxScaleWidTDLU=$voronoi_maxScaleWidTDLU
# minimum strength
minStrTDLU=$voronoi_minStrTDLU
# maximum strength
maxStrTDLU=$voronoi_maxStrTDLU
# maximum deflection (fraction of pi)
fatInFatDeflectMax=$voronoi_fatInFatDeflectMax
# minimum length scale
minScaleLenFatInFat=$voronoi_minScaleLenFatInFat
# maximum length scale
maxScaleLenFatInFat=$voronoi_maxScaleLenFatInFat
# minimum width scale
minScaleWidFatInFat=$voronoi_minScaleWidFatInFat
# maximum width scale
maxScaleWidFatInFat=$voronoi_maxScaleWidFatInFat
# minimum strength
minStrFatInFat=$voronoi_minStrFatInFat
# maximum strength
maxStrFatInFat=$voronoi_maxStrFatInFat
# maximum deflection (fraction of pi)
fatInGlandDeflectMax=$voronoi_fatInGlandDeflectMax
# minimum length scale
minScaleLenFatInGland=$voronoi_minScaleLenFatInGland
# maximum length scale
maxScaleLenFatInGland=$voronoi_maxScaleLenFatInGland
# minimum width scale
minScaleWidFatInGland=$voronoi_minScaleWidFatInGland
# maximum width scale
maxScaleWidFatInGland=$voronoi_maxScaleWidFatInGland
# minimum strength
minStrFatInGland=$voronoi_minStrFatInGland
# maximum strength
maxStrFatInGland=$voronoi_maxStrFatInGland
# maximum deflection (fraction of pi)
glandInGlandDeflectMax=$voronoi_glandInGlandDeflectMax
# minimum length scale
minScaleLenGlandInGland=$voronoi_minScaleLenGlandInGland
# maximum length scale
maxScaleLenGlandInGland=$voronoi_maxScaleLenGlandInGland
# minimum width scale
minScaleWidGlandInGland=$voronoi_minScaleWidGlandInGland
# maximum width scale
maxScaleWidGlandInGland=$voronoi_maxScaleWidGlandInGland
# minimum strength
minStrGlandInGland=$voronoi_minStrGlandInGland
# maximum strength
maxStrGlandInGland=$voronoi_maxStrGlandInGland
# check seeds in radius (mm) 
seedRadius=$voronoi_seedRadius

#####################
# fat variables
#####################
[fat]
# min lobule axis length (mm)
# was 20
minLobuleAxis=$fat_minLobuleAxis
# max lobule axis length (mm)
# was 35
maxLobuleAxis=$fat_maxLobuleAxis
# axial ratio min
minAxialRatio=$fat_minAxialRatio
# axial ratio max
maxAxialRatio=$fat_maxAxialRatio
# minimum ligament separation between lobules
minLobuleGap=$fat_minLobuleGap
# maximum of absolute value of Fourier coefficient as fraction of main radius
maxCoeffStr=$fat_maxCoeffStr
# minimum of absolute value of Fourier coefficient as fraction of main radius
minCoeffStr=$fat_minCoeffStr
# maximum number of trial lobules
# was 401
maxLobuleTry=$fat_maxLobuleTry


#####################
# ligament variables
#####################
[lig]
thickness=$ligament_thickness
targetFrac=$ligament_targetFrac
maxTry=$ligament_maxTry
minAxis=$ligament_minAxis
maxAxis=$ligament_maxAxis
minAxialRatio=$ligament_minAxialRatio
maxAxialRatio=$ligament_maxAxialRatio
maxPerturb=$ligament_maxPerturb
maxDeflect=$ligament_maxDeflect
scale=$ligament_scale
lacunarity=$ligament_lacunarity
persistence=$ligament_persistence
numOctaves=$ligament_numOctaves

#####################
# duct tree variables
#####################
[ductTree]
# target number of branches (uint)
maxBranch=$duct_tree_maxBranch
# maximum generation (uint)
maxGen=$duct_tree_maxGen
# initial radius of tree (mm)
initRad=$duct_tree_initRad
# base Length of root duct at nipple (mm)
baseLength=$duct_tree_baseLength
# number of voxels for tree density tracking (uint)
nFillX=$duct_tree_nFillX
nFillY=$duct_tree_nFillY
nFillZ=$duct_tree_nFillZ

#####################
# duct branch variables
#####################
[ductBr]
# minimum branch radius to have children (mm)
childMinRad=$duct_branch_childMinRad
# minimum starting radius as a fraction of parent end radius
minRadFrac=$duct_branch_minRadFrac
# maximum starting radius as a fraction of parent end radius
maxRadFrac=$duct_branch_maxRadFrac
# length reduction as fraction of parent length
lenShrink=$duct_branch_lenShrink
# maximum jitter in branch length (fraction)
lenRange=$duct_branch_lenRange
# aximuthal angle noise (radians)
rotateJitter=$duct_branch_rotateJitter

#####################
# duct segment variables
#####################
[ductSeg]
# radius distribution shape parameters
radiusBetaA=$duct_segment_radiusBetaA
radiusBetaB=$duct_segment_radiusBetaB
# fraction of branch length per segment
segFrac=$duct_segment_segFrac
# maximum radius of curvature (mm)
maxCurvRad=$duct_segment_maxCurvRad
# maximum length of segment based on
# curvature (fraction of pi radians)
maxCurvFrac=$duct_segment_maxCurvFrac
# min and max end radius as fraction of start radius
minEndRad=$duct_segment_minEndRad
maxEndRad=$duct_segment_maxEndRad
# cost function preferential angle weighting
angleWt=$duct_segment_angleWt
# cost function density weighting
densityWt=$duct_segment_densityWt
# number of trial segments to generate (uint)
numTry=$duct_segment_numTry
# maximum number of segments to generate before
# giving up and reducing length (uint)
maxTry=$duct_segment_maxTry
# total number of segment tries before completely giving up
absMaxTry=$duct_segment_absMaxTry
# step size for checking segment is valid (mm)
roiStep=$duct_segment_roiStep


#####################
# vessel tree variables
#####################
[vesselTree]
# target number of branches (uint)
maxBranch=$vessel_tree_maxBranch
# maximum generation (uint)
maxGen=$vessel_tree_maxGen
# initial radius of tree (mm)
initRad=$vessel_tree_initRad
# base length of root vessel (mm)
baseLength=$vessel_tree_baseLength
# number of voxels for tree density tracking (uint)
nFillX=$vessel_tree_nFillX
nFillY=$vessel_tree_nFillY
nFillZ=$vessel_tree_nFillZ

#####################
# vessel branch variables
#####################
[vesselBr]
# minimum branch radius to have children (mm)
childMinRad=$vessel_branch_childMinRad
# minimum starting radius as a fraction of parent end radius 
minRadFrac=$vessel_branch_minRadFrac
# maximum starting radius as a fraction of parent end radius 
maxRadFrac=$vessel_branch_maxRadFrac
# length reduction as fraction of parent length
lenShrink=$vessel_branch_lenShrink
# maximum jitter in branch length (fraction)
lenRange=$vessel_branch_lenRange
# aximuthal angle noise (radians)
rotateJitter=$vessel_branch_rotateJitter

#####################
# vessel segment variables
#####################
[vesselSeg]
# radius distribution shape parameters
radiusBetaA=$vessel_segment_radiusBetaA
radiusBetaB=$vessel_segment_radiusBetaB
# fraction of branch length to segment
segFrac=$vessel_segment_segFrac
# maximum radius of curvature (mm)
maxCurvRad=$vessel_segment_maxCurvRad
# maximum length of segment based on 
# curvature (fraction of pi radians)
maxCurvFrac=$vessel_segment_maxCurvFrac
# min and max end radius as fraction of start radius
minEndRad=$vessel_segment_minEndRad
maxEndRad=$vessel_segment_maxEndRad
# cost function preferential angle weighting
angleWt=$vessel_segment_angleWt
# cost function density weighting
densityWt=$vessel_segment_densityWt
# cost function direction weighting
dirWt=$vessel_segment_dirWt
# number of trial segments to generate (uint)
numTry=$vessel_segment_numTry
# maximum number of segments to generate before
# giving up and reducing length (uint)
maxTry=$vessel_segment_maxTry
# total number of segment tries before completely giving up
absMaxTry=$vessel_segment_absMaxTry
# step size for checking segment is valid (mm)
roiStep=$vessel_segment_roiStep