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
ures=$ures
# v resolution of base shape
vres=$vres
# minimum point separation (mm)
pointSep=$pointSep
# back ring thickness (mm)
ringWidth=$ringWidth
# back ring step size (mm)
ringSep=$ringSep
# angle to preserve while smoothing (degrees)
featureAngle=$featureAngle
# fraction of triangles to decimate
targetReduction=$targetReduction
# bottom scale was 1.3
a1b=$a1b
# top scale was 1.3
a1t=$a1t
# left scale
a2l=$a2l
# right scale
a2r=$a2r
# outward scale was 1.75
a3=$a3
# u quadric exponent
eps1=$eps1
# v quadric exponent
eps2=$eps2
# do ptosis deformation (boolean)
doPtosis=$doPtosis
ptosisB0=$ptosisB0
ptosisB1=$ptosisB1
# do turn deformation (boolean)
doTurn=$doTurn
turnC0=$turnC0
turnC1=$turnC1
# do top shape deformation (boolean)
doTopShape=$doTopShape
topShapeS0=$topShapeS0
topShapeS1=$topShapeS1
topShapeT0=$topShapeT0
topShapeT1=$topShapeT1
# do flatten size deformation (boolean)
doFlattenSide=$doFlattenSide
flattenSideG0=$flattenSideG0
flattenSideG1=$flattenSideG1
# do turn top deformation (boolean)
doTurnTop=$doTurnTop
turnTopH0=$turnTopH0
turnTopH1=$turnTopH1

#####################
# breast compartment
#####################
[compartments]
# number of breast compartments
num=$num
# distance along nipple line of compartment seed base (mm)
seedBaseDist=$seedBaseDist
# fraction of phantom in nipple direction forced to be fat
backFatBufferFrac=$backFatBufferFrac
# number of backplane seed points
numBackSeeds=$numBackSeeds
# maximum seed jitter (fraction of subtended angle)
angularJitter=$angularJitter
# maximum seed jitter in nipple direction (mm)
zJitter=$zJitter
# maximum radial distance from base seed as a fraction of distance to breast surface
maxFracRadialDist=$maxFracRadialDist
# minimum radial distance from base seed as a fraction of distance to breast surface
minFracRadialDist=$minFracRadialDist
# minimum scale in nipple direction
minScaleNippleDir=$minScaleNippleDir
# maximum scale in nipple direction
maxScaleNippleDir=$maxScaleNippleDir
# minimum scale in non-nipple direction
minScale=$minScale
# maximum scale in non-nipple direction
maxScale=$maxScale
# minimum gland strength
minGlandStrength=$minGlandStrength
# maximum gland strength
maxGlandStrength=$maxGlandStrength
# maximum compartment deflection angle from pointing towards nipple (fraction of pi)
maxDeflect=$maxDeflect
# minimum scale skin seeds in nipple direction
minSkinScaleNippleDir=$minSkinScaleNippleDir
# maximum scale skin seeds in nipple direction
maxSkinScaleNippleDir=$maxSkinScaleNippleDir
# minimum scale skin in non-nipple direction
minSkinScale=$minSkinScale
# maximum scale skin in non-nipple direction
maxSkinScale=$maxSkinScale
# skin strength
skinStrength=$skinStrength
# back scale
backScale=$backScale
# back strength
backStrength=$backStrength
# nipple scale
nippleScale=$nippleScale
# nipple strength
nippleStrength=$nippleStrength
# check seeds within radius (mm)
voronSeedRadius=$voronSeedRadius

#####################
# TDLU variables
#####################
[TDLU]
# maximum TDLU length
maxLength=$maxLength
# minimum TDLU length
minLength=$minLength
# maximum TDLU width
maxWidth=$maxWidth
# minimum TDLU width
minWidth=$minWidth

#####################
# Perlin noise variables
#####################
[perlin]
# maximum fraction of radius deviation 
maxDeviation=$maxDeviation
# starting frequency
frequency=$frequency
# octave frequency multiplier
lacunarity=$lacunarity
# octave signal decay
persistence=$persistence
# number of frequency octaves
numOctaves=$numOctaves
# x direction noise generation seed
xNoiseGen=$xNoiseGen
# y direction noise generation seed
yNoiseGen=$yNoiseGen
# z direction noise generation seed
zNoiseGen=$zNoiseGen
# seed noise generation
seedNoiseGen=$seedNoiseGen
# shift noise generation seed
shiftNoiseGen=$shiftNoiseGen

#####################
# Compartment boundary noise
#####################
[boundary]
# maximum fraction of distance deviation 
maxDeviation=$maxDeviation
# starting frequency
frequency=$frequency
# octave frequency multiplier
lacunarity=$lacunarity
# octave signal decay
persistence=$persistence

#####################
# Lobule boundary perturbation noise
#####################
[perturb]
# maximum fraction of distance deviation 
maxDeviation=$maxDeviation
# starting frequency
frequency=$frequency
# octave frequency multiplier
lacunarity=$lacunarity
# octave signal decay
persistence=$persistence

#####################
# Lobule glandular buffer noise
#####################
[buffer]
# maximum fraction of distance deviation 
maxDeviation=$maxDeviation
# starting frequency
frequency=$frequency
# octave frequency multiplier
lacunarity=$lacunarity
# octave signal decay
persistence=$persistence

#####################
# Voronoi segmentation variables
#####################
[voronoi]
# fat voronoi seed density (mm^-3)
fatInFatSeedDensity=$fatInFatSeedDensity
# fat voronoi seed in glandular tissue density (mm^-3)
fatInGlandSeedDensity=$fatInGlandSeedDensity
# glandular voronoi seed density (mm^-3)
glandInGlandSeedDensity=$glandInGlandSeedDensity
# maximum deflection (fraction of pi)
TDLUDeflectMax=$TDLUDeflectMax
# minimum length scale
minScaleLenTDLU=$minScaleLenTDLU
# maximum length scale
maxScaleLenTDLU=$maxScaleLenTDLU
# minimum width scale
minScaleWidTDLU=$minScaleWidTDLU
# maximum width scale
maxScaleWidTDLU=$maxScaleWidTDLU
# minimum strength
minStrTDLU=$minStrTDLU
# maximum strength
maxStrTDLU=$maxStrTDLU
# maximum deflection (fraction of pi)
fatInFatDeflectMax=$fatInFatDeflectMax
# minimum length scale
minScaleLenFatInFat=$minScaleLenFatInFat
# maximum length scale
maxScaleLenFatInFat=$maxScaleLenFatInFat
# minimum width scale
minScaleWidFatInFat=$minScaleWidFatInFat
# maximum width scale
maxScaleWidFatInFat=$maxScaleWidFatInFat
# minimum strength
minStrFatInFat=$minStrFatInFat
# maximum strength
maxStrFatInFat=$maxStrFatInFat
# maximum deflection (fraction of pi)
fatInGlandDeflectMax=$fatInGlandDeflectMax
# minimum length scale
minScaleLenFatInGland=$minScaleLenFatInGland
# maximum length scale
maxScaleLenFatInGland=$maxScaleLenFatInGland
# minimum width scale
minScaleWidFatInGland=$minScaleWidFatInGland
# maximum width scale
maxScaleWidFatInGland=$maxScaleWidFatInGland
# minimum strength
minStrFatInGland=$minStrFatInGland
# maximum strength
maxStrFatInGland=$maxStrFatInGland
# maximum deflection (fraction of pi)
glandInGlandDeflectMax=$glandInGlandDeflectMax
# minimum length scale
minScaleLenGlandInGland=$minScaleLenGlandInGland
# maximum length scale
maxScaleLenGlandInGland=$maxScaleLenGlandInGland
# minimum width scale
minScaleWidGlandInGland=$minScaleWidGlandInGland
# maximum width scale
maxScaleWidGlandInGland=$maxScaleWidGlandInGland
# minimum strength
minStrGlandInGland=$minStrGlandInGland
# maximum strength
maxStrGlandInGland=$maxStrGlandInGland
# check seeds in radius (mm) 
seedRadius=$seedRadius

#####################
# fat variables
#####################
[fat]
# min lobule axis length (mm)
# was 20
minLobuleAxis=$minLobuleAxis
# max lobule axis length (mm)
# was 35
maxLobuleAxis=$maxLobuleAxis
# axial ratio min
minAxialRatio=$minAxialRatio
# axial ratio max
maxAxialRatio=$maxAxialRatio
# minimum ligament separation between lobules
minLobuleGap=$minLobuleGap
# maximum of absolute value of Fourier coefficient as fraction of main radius
maxCoeffStr=$maxCoeffStr
# minimum of absolute value of Fourier coefficient as fraction of main radius
minCoeffStr=$minCoeffStr
# maximum number of trial lobules
# was 401
maxLobuleTry=$maxLobuleTry


#####################
# ligament variables
#####################
[lig]
thickness=$thickness
targetFrac=$targetFrac
maxTry=$maxTry
minAxis=$minAxis
maxAxis=$maxAxis
minAxialRatio=$minAxialRatio
maxAxialRatio=$maxAxialRatio
maxPerturb=$maxPerturb
maxDeflect=$maxDeflect
scale=$scale
lacunarity=$lacunarity
persistence=$persistence
numOctaves=$numOctaves

#####################
# duct tree variables
#####################
[ductTree]
# target number of branches (uint)
maxBranch=$maxBranch
# maximum generation (uint)
maxGen=$maxGen
# initial radius of tree (mm)
initRad=$initRad
# base Length of root duct at nipple (mm)
baseLength=$baseLength
# number of voxels for tree density tracking (uint)
nFillX=$nFillX
nFillY=$nFillY
nFillZ=$nFillZ

#####################
# duct branch variables
#####################
[ductBr]
# minimum branch radius to have children (mm)
childMinRad=$childMinRad
# minimum starting radius as a fraction of parent end radius
minRadFrac=$minRadFrac
# maximum starting radius as a fraction of parent end radius
maxRadFrac=$maxRadFrac
# length reduction as fraction of parent length
lenShrink=$lenShrink
# maximum jitter in branch length (fraction)
lenRange=$lenRange
# aximuthal angle noise (radians)
rotateJitter=$rotateJitter

#####################
# duct segment variables
#####################
[ductSeg]
# radius distribution shape parameters
radiusBetaA=$radiusBetaA
radiusBetaB=$radiusBetaB
# fraction of branch length per segment
segFrac=$segFrac
# maximum radius of curvature (mm)
maxCurvRad=$maxCurvRad
# maximum length of segment based on
# curvature (fraction of pi radians)
maxCurvFrac=$maxCurvFrac
# min and max end radius as fraction of start radius
minEndRad=$minEndRad
maxEndRad=$maxEndRad
# cost function preferential angle weighting
angleWt=$angleWt
# cost function density weighting
densityWt=$densityWt
# number of trial segments to generate (uint)
numTry=$numTry
# maximum number of segments to generate before
# giving up and reducing length (uint)
maxTry=$maxTry
# total number of segment tries before completely giving up
absMaxTry=$absMaxTry
# step size for checking segment is valid (mm)
roiStep=$roiStep


#####################
# vessel tree variables
#####################
[vesselTree]
# target number of branches (uint)
maxBranch=$maxBranch
# maximum generation (uint)
maxGen=$maxGen
# initial radius of tree (mm)
initRad=$initRad
# base length of root vessel (mm)
baseLength=$baseLength
# number of voxels for tree density tracking (uint)
nFillX=$nFillX
nFillY=$nFillY
nFillZ=$nFillZ

#####################
# vessel branch variables
#####################
[vesselBr]
# minimum branch radius to have children (mm)
childMinRad=$childMinRad
# minimum starting radius as a fraction of parent end radius 
minRadFrac=$minRadFrac
# maximum starting radius as a fraction of parent end radius 
maxRadFrac=$maxRadFrac
# length reduction as fraction of parent length
lenShrink=$lenShrink
# maximum jitter in branch length (fraction)
lenRange=$lenRange
# aximuthal angle noise (radians)
rotateJitter=$rotateJitter

#####################
# vessel segment variables
#####################
[vesselSeg]
# radius distribution shape parameters
radiusBetaA=$radiusBetaA
radiusBetaB=$radiusBetaB
# fraction of branch length to segment
segFrac=$segFrac
# maximum radius of curvature (mm)
maxCurvRad=$maxCurvRad
# maximum length of segment based on 
# curvature (fraction of pi radians)
maxCurvFrac=$maxCurvFrac
# min and max end radius as fraction of start radius
minEndRad=$minEndRad
maxEndRad=$maxEndRad
# cost function preferential angle weighting
angleWt=$angleWt
# cost function density weighting
densityWt=$densityWt
# cost function direction weighting
dirWt =$dirWt 
# number of trial segments to generate (uint)
numTry=$numTry
# maximum number of segments to generate before
# giving up and reducing length (uint)
maxTry=$maxTry
# total number of segment tries before completely giving up
absMaxTry=$absMaxTry
# step size for checking segment is valid (mm)
roiStep=$roiStep