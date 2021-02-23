/*
 * breastPhantom.cxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

// create volumetric breast

#include "breastPhantom.hxx"

namespace po = boost::program_options;

unsigned int ductTree::num = 0;
unsigned int arteryTree::num = 0;
unsigned int veinTree::num = 0;

int main(int argc, char* argv[]){

  double pi = vtkMath::Pi();
	
  // timing
  //time_t currentTime, previousTime;
  //time(&previousTime);
  //time_t startTime, endTime;
  //time(&startTime);

  // configuration file variables
  po::options_description baseOpt("base options");
  baseOpt.add_options()
    ("base.outputDir",po::value<std::string>()->default_value("."),"output directory")
    ("base.imgRes",po::value<double>()->default_value(0.25),"phantom resolution (mm)")
    ("base.skinThick",po::value<double>()->default_value(0.75),"skin thickness (mm)")
    ("base.nippleLen",po::value<double>()->default_value(4.0),"nipple length (mm)")
    ("base.nippleRad",po::value<double>()->default_value(4.0),"nipple radius (mm)")
    ("base.areolaRad",po::value<double>()->default_value(8.0),"areola radius (mm)")
    ("base.leftBreast",po::value<bool>()->default_value(true),"left side breast (boolean)")
    ("base.targetFatFrac",po::value<double>()->default_value(0.75),"desired fraction of breast to be fat")
    ("base.seed",po::value<unsigned int>(),"random number generator seed")
    ;

  po::options_description shapeOpt("breast shape options");
  shapeOpt.add_options()
    ("shape.ures",po::value<double>()->default_value(0.02),"u resolution of base shape")
    ("shape.vres",po::value<double>()->default_value(0.02),"v resolution of base shape")
    ("shape.pointSep",po::value<double>()->default_value(0.01),"minimum point separation (mm)")
    ("shape.ringWidth",po::value<double>()->default_value(10.0),"thickness of back ring (mm)")
    ("shape.ringSep",po::value<double>()->default_value(0.5),"ring node step size (mm)")
    ("shape.featureAngle",po::value<double>()->default_value(20.0),"angle to preserve while smoothing (degrees)")
    ("shape.targetReduction",po::value<double>()->default_value(0.05),"fraction of triangles to decimate")
    ("shape.a1b",po::value<double>()->default_value(1.0),"bottom size")
    ("shape.a1t",po::value<double>()->default_value(1.0),"top size")
    ("shape.a2l",po::value<double>()->default_value(1.0),"left size")
    ("shape.a2r",po::value<double>()->default_value(1.0),"right size")
    ("shape.a3",po::value<double>()->default_value(1.0),"outward size")
    ("shape.eps1",po::value<double>()->default_value(1.0),"u quadric exponent")
    ("shape.eps2",po::value<double>()->default_value(1.0),"v quadric exponent")
    ("shape.doPtosis",po::value<bool>()->default_value(false),"ptosis deformation (boolean)")
    ("shape.ptosisB0",po::value<double>()->default_value(0.4),"ptosis b0")
    ("shape.ptosisB1",po::value<double>()->default_value(0.3),"ptosis b1")
    ("shape.doTurn",po::value<bool>()->default_value(false),"turn deformation (boolean)")
    ("shape.turnC0",po::value<double>()->default_value(-0.4),"turn c0")
    ("shape.turnC1",po::value<double>()->default_value(0.2),"turn c1")
    ("shape.doTopShape",po::value<bool>()->default_value(false),"top shape deformation (boolean)")
    ("shape.topShapeS0",po::value<double>()->default_value(1.54),"top shape s0")
    ("shape.topShapeS1",po::value<double>()->default_value(1.3),"top shape s1")
    ("shape.topShapeT0",po::value<double>()->default_value(-3.0),"top shape t0")
    ("shape.topShapeT1",po::value<double>()->default_value(-1.0),"top shape t1")
    ("shape.doFlattenSide",po::value<bool>()->default_value(false),"flatten side deformation (boolean)")
    ("shape.flattenSideG0",po::value<double>()->default_value(-0.4),"flatten side g0")
    ("shape.flattenSideG1",po::value<double>()->default_value(0.2),"flatten side g1")
    ("shape.doTurnTop",po::value<bool>()->default_value(false),"turn top deformation (boolean)")
    ("shape.turnTopH0",po::value<double>()->default_value(-2.7),"turn top h0")
    ("shape.turnTopH1",po::value<double>()->default_value(-3.7),"turn top h1")
    ;

  po::options_description compartOpt("breast compartment options");
  compartOpt.add_options()
    ("compartments.num",po::value<int>()->default_value(8),"number of breast compartments")
    ("compartments.seedBaseDist",po::value<double>()->default_value(12.5),"distance along nipple line of compartment seed base (mm)")
    ("compartments.backFatBufferFrac",po::value<double>()->default_value(0.2),"fraction of phantom in nipple direction forced to be fat")
    ("compartments.numBackSeeds",po::value<int>()->default_value(100),"number of backplane seed points")
    ("compartments.angularJitter",po::value<double>()->default_value(0.125),"maximum seed jitter (fraction of subtended angle)")
    ("compartments.zJitter",po::value<double>()->default_value(5.0),"maximum seed jitter in nipple direction (mm)")
    ("compartments.maxFracRadialDist",po::value<double>()->default_value(0.5),"maximum radial distance from base seed as a fraction of distance to breast surface")
    ("compartments.minFracRadialDist",po::value<double>()->default_value(0.25),"minimum radial distance from base seed as a fraction of distance to breast surface")
    ("compartments.minScaleNippleDir",po::value<double>()->default_value(0.01),"minimum scale in nipple direction")
    ("compartments.maxScaleNippleDir",po::value<double>()->default_value(0.01),"maximum scale in nipple direction")
    ("compartments.minScale",po::value<double>()->default_value(50.0),"minimum scale in non-nipple direction")
    ("compartments.maxScale",po::value<double>()->default_value(60.0),"maximum scale in non-nipple direction")
    ("compartments.minGlandStrength",po::value<double>()->default_value(20.0),"minimum gland strength")
    ("compartments.maxGlandStrength",po::value<double>()->default_value(20.0),"maximum gland strength")
    ("compartments.maxDeflect",po::value<double>()->default_value(0.01),"maximum compartment deflection angle from pointing towards nipple (fraction of pi)")
    ("compartments.minSkinScaleNippleDir",po::value<double>()->default_value(0.01),"minimum scale skin seeds in nipple direction")
    ("compartments.maxSkinScaleNippleDir",po::value<double>()->default_value(0.01),"maximum scale skin seeds in nipple direction")
    ("compartments.minSkinScale",po::value<double>()->default_value(70.0),"minimum scale skin in non-nipple direction")
    ("compartments.maxSkinScale",po::value<double>()->default_value(70.0),"maximum scale skin in non-nipple direction")
    ("compartments.skinStrength",po::value<double>()->default_value(20.0),"skin strength")
    ("compartments.backScale",po::value<double>()->default_value(50.0),"back scale")
    ("compartments.backStrength",po::value<double>()->default_value(50.0),"back strength")
    ("compartments.nippleScale",po::value<double>()->default_value(50.0),"nipple scale")
    ("compartments.nippleStrength",po::value<double>()->default_value(50.0),"nipple strength")
    ("compartments.voronSeedRadius",po::value<double>()->default_value(40.0),"check all seeds in radius (mm)")
    ;
  
  po::options_description TDLUOpt("TDLU options");
  TDLUOpt.add_options()
    ("TDLU.maxLength",po::value<double>()->default_value(4.0),"maximum TDLU length")
    ("TDLU.minLength",po::value<double>()->default_value(2.0),"minimum TDLU length")
    ("TDLU.maxWidth",po::value<double>()->default_value(2.0),"maximum TDLU width")
    ("TDLU.minWidth",po::value<double>()->default_value(1.0),"minimum TDLU width")
    ;
  
  po::options_description perlinOpt("Perlin noise options");
  perlinOpt.add_options()
    ("perlin.maxDeviation",po::value<double>()->default_value(0.1),"maximum perlin perturbation fraction of radius")
    ("perlin.frequency",po::value<double>()->default_value(0.5),"starting frequency")
    ("perlin.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
    ("perlin.persistence",po::value<double>()->default_value(0.2),"octave signal decay")
    ("perlin.numOctaves",po::value<int>()->default_value(6),"number of frequency octaves")
    ("perlin.xNoiseGen",po::value<int>()->default_value(683),"x direction noise generation seed")
    ("perlin.yNoiseGen",po::value<int>()->default_value(4933),"y direction noise generation seed")
    ("perlin.zNoiseGen",po::value<int>()->default_value(23),"z direction noise generation seed")
    ("perlin.seedNoiseGen",po::value<int>()->default_value(3095),"seed noise generation")
    ("perlin.shiftNoiseGen",po::value<int>()->default_value(11),"shift noise generation seed")
    ;

  po::options_description boundaryOpt("Boundary noise options");
  boundaryOpt.add_options()
    ("boundary.maxDeviation",po::value<double>()->default_value(0.25),"maximum perlin perturbation fraction of radius")
    ("boundary.frequency",po::value<double>()->default_value(0.2),"starting frequency")
    ("boundary.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
    ("boundary.persistence",po::value<double>()->default_value(0.5),"octave signal decay")
    ;

  po::options_description perturbOpt("Lobule perturbation noise options");
  perturbOpt.add_options()
    ("perturb.maxDeviation",po::value<double>()->default_value(0.25),"maximum perlin perturbation fraction of radius")
    ("perturb.frequency",po::value<double>()->default_value(0.2),"starting frequency")
    ("perturb.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
    ("perturb.persistence",po::value<double>()->default_value(0.5),"octave signal decay")
    ;
	
  po::options_description bufferOpt("Lobule glandular buffer options");
  bufferOpt.add_options()
    ("buffer.maxDeviation",po::value<double>()->default_value(0.1),"maximum perlin buffer fraction of radius")
    ("buffer.frequency",po::value<double>()->default_value(0.2),"starting frequency")
    ("buffer.lacunarity",po::value<double>()->default_value(1.5),"octave frequency multiplier")
    ("buffer.persistence",po::value<double>()->default_value(0.5),"octave signal decay")
    ;
	
  po::options_description voronOpt("Voronoi options");
  voronOpt.add_options()
    ("voronoi.fatInFatSeedDensity",po::value<double>()->default_value(0.0),"fat voronoi seed density (mm^-3)")
    ("voronoi.fatInGlandSeedDensity",po::value<double>()->default_value(0.0),"fat voronoi seed in glandular tissue density (mm^-3)")
    ("voronoi.glandInGlandSeedDensity",po::value<double>()->default_value(0.0),"glandular voronoi seed density (mm^-3)")
    ("voronoi.TDLUDeflectMax",po::value<double>()->default_value(0.0),"maximum deflection (fraction of pi)")
    ("voronoi.minScaleLenTDLU",po::value<double>()->default_value(0.0),"minimum length scale")
    ("voronoi.maxScaleLenTDLU",po::value<double>()->default_value(0.0),"maximum length scale")
    ("voronoi.minScaleWidTDLU",po::value<double>()->default_value(0.0),"minimum width scale")
    ("voronoi.maxScaleWidTDLU",po::value<double>()->default_value(0.0),"maximum width scale")
    ("voronoi.minStrTDLU",po::value<double>()->default_value(0.0),"minimum strength")
    ("voronoi.maxStrTDLU",po::value<double>()->default_value(0.0),"maximum strength")
    ("voronoi.fatInFatDeflectMax",po::value<double>()->default_value(0.15),"maximum deflection (fraction of pi)")
    ("voronoi.minScaleLenFatInFat",po::value<double>()->default_value(0.0),"minimum length scale")
    ("voronoi.maxScaleLenFatInFat",po::value<double>()->default_value(0.0),"maximum length scale")
    ("voronoi.minScaleWidFatInFat",po::value<double>()->default_value(0.0),"minimum width scale")
    ("voronoi.maxScaleWidFatInFat",po::value<double>()->default_value(0.0),"maximum width scale")
    ("voronoi.minStrFatInFat",po::value<double>()->default_value(0.0),"minimum strength")
    ("voronoi.maxStrFatInFat",po::value<double>()->default_value(0.0),"maximum strength")
    ("voronoi.fatInGlandDeflectMax",po::value<double>()->default_value(0.0),"maximum deflection (fraction of pi)")
    ("voronoi.minScaleLenFatInGland",po::value<double>()->default_value(0.0),"minimum length scale")
    ("voronoi.maxScaleLenFatInGland",po::value<double>()->default_value(0.0),"maximum length scale")
    ("voronoi.minScaleWidFatInGland",po::value<double>()->default_value(0.0),"minimum width scale")
    ("voronoi.maxScaleWidFatInGland",po::value<double>()->default_value(0.0),"maximum width scale")
    ("voronoi.minStrFatInGland",po::value<double>()->default_value(0.0),"minimum strength")
    ("voronoi.maxStrFatInGland",po::value<double>()->default_value(0.0),"maximum strength")
    ("voronoi.glandInGlandDeflectMax",po::value<double>()->default_value(0.0),"maximum deflection (fraction of pi)")
    ("voronoi.minScaleLenGlandInGland",po::value<double>()->default_value(0.0),"minimum length scale")
    ("voronoi.maxScaleLenGlandInGland",po::value<double>()->default_value(0.0),"maximum length scale")
    ("voronoi.minScaleWidGlandInGland",po::value<double>()->default_value(0.0),"minimum width scale")
    ("voronoi.maxScaleWidGlandInGland",po::value<double>()->default_value(0.0),"maximum width scale")
    ("voronoi.minStrGlandInGland",po::value<double>()->default_value(0.0),"minimum strength")
    ("voronoi.maxStrGlandInGland",po::value<double>()->default_value(0.0),"maximum strength")
    ("voronoi.seedRadius",po::value<double>()->default_value(10),"check seeds in radius (mm)")
    ;

  po::options_description fatOpt("Fat lobule options");
  fatOpt.add_options()
    ("fat.minLobuleAxis",po::value<double>()->default_value(7.0),"min lobule axis length (mm)")
    ("fat.maxLobuleAxis",po::value<double>()->default_value(12.0),"max lobule axis length (mm)")
    ("fat.minAxialRatio",po::value<double>()->default_value(0.33),"min axial ratio")
    ("fat.maxAxialRatio",po::value<double>()->default_value(0.53),"max axial ratio")
    ("fat.minLobuleGap",po::value<double>()->default_value(0.2),"minimum ligament separation between lobules")
    ("fat.maxCoeffStr",po::value<double>()->default_value(0.08),"maximum of absolute value of Fourier coefficient as fraction of main radius")
    ("fat.minCoeffStr",po::value<double>()->default_value(0.02),"minimum of absolute value of Fourier coefficient as fraction of main radius")
    ("fat.maxLobuleTry",po::value<int>()->default_value(600),"maximum number of trial lobules")
    ;

  po::options_description ligOpt("Ligament options");
  ligOpt.add_options()
    ("lig.thickness",po::value<double>()->default_value(0.12),"ligament thickness (mm)")
    ("lig.targetFrac",po::value<double>()->default_value(0.85),"ligamented volume stopping fraction")
    ("lig.maxTry",po::value<int>()->default_value(15000),"maximum number of ligaments")
    ("lig.minAxis",po::value<double>()->default_value(15.0),"min lobule axis length (mm)")
    ("lig.maxAxis",po::value<double>()->default_value(20.0),"max lobule axis length (mm)")
    ("lig.minAxialRatio",po::value<double>()->default_value(0.3),"min axial ratio")
    ("lig.maxAxialRatio",po::value<double>()->default_value(0.4),"max axial ratio")
    ("lig.maxPerturb",po::value<double>()->default_value(0.1),"maximum perlin perturbation fraction of radius")
    ("lig.maxDeflect",po::value<double>()->default_value(0.12),"maximum deflection from nipple direction (fraction of pi)")
    ("lig.scale",po::value<double>()->default_value(0.007),"perlin frequency scaling")
    ("lig.lacunarity",po::value<double>()->default_value(1.8),"octave frequency multiplier")
    ("lig.persistence",po::value<double>()->default_value(0.6),"octave signal decay")
    ("lig.numOctaves",po::value<int>()->default_value(6),"number of frequency octaves")
    ;

  po::options_description ductTreeOpt("Duct tree options");
  ductTreeOpt.add_options()
    ("ductTree.maxBranch",po::value<uint>()->default_value(100),"Maximum number of branches")
    ("ductTree.maxGen",po::value<uint>()->default_value(15),"Maximum generation")
    ("ductTree.baseLength",po::value<double>()->default_value(7.6),"main branch length (mm)")
    ("ductTree.initRad",po::value<double>()->default_value(2.0),"tree start radius")
    ("ductTree.nFillX",po::value<uint>()->default_value(100),"number x voxels for density map")
    ("ductTree.nFillY",po::value<uint>()->default_value(100),"number y voxels for density map")
    ("ductTree.nFillZ",po::value<uint>()->default_value(100),"number z voxels for density map")
    ;

  po::options_description ductBrOpt("Duct Branch options");
  ductBrOpt.add_options()
    ("ductBr.childMinRad",po::value<double>()->default_value(0.25),"min radius to have children")
    ("ductBr.minRadFrac",po::value<double>()->default_value(0.65),"min start radius as fraction of parent radius")
    ("ductBr.maxRadFrac",po::value<double>()->default_value(0.99),"max start radius as fraction of parent radius")
    ("ductBr.lenShrink",po::value<double>()->default_value(0.5),"duct length relative to parent")
    ("ductBr.lenRange",po::value<double>()->default_value(0.1),"duct length fraction variability")
    ("ductBr.rotateJitter",po::value<double>()->default_value(0.1),"jitter in azimuthal angle of second child")
    ;
  
  po::options_description ductSegOpt("Duct Segment options");
  ductSegOpt.add_options()
    ("ductSeg.radiusBetaA",po::value<double>()->default_value(2.0),"radius distribution shape parameter A")
    ("ductSeg.radiusBetaB",po::value<double>()->default_value(2.0),"radius distribution shape parameter B")
    ("ductSeg.maxCurvRad",po::value<double>()->default_value(20.0),"max radius of curvature")
    ("ductSeg.maxCurvFrac",po::value<double>()->default_value(0.33),"max fraction of circle")
    ("ductSeg.minEndRad",po::value<double>()->default_value(0.85),"min end radius as fraction of start radius")
    ("ductSeg.maxEndRad",po::value<double>()->default_value(1.05),"max end radius as fraction of start radius")
    ("ductSeg.angleWt",po::value<double>()->default_value(1.0),"cost function preferential angle weighting")
    ("ductSeg.densityWt",po::value<double>()->default_value(5e-5),"cost function density weighting")
    ("ductSeg.numTry",po::value<uint>()->default_value(10),"number of trial segments")
    ("ductSeg.maxTry",po::value<uint>()->default_value(100),"max number of trial segments before reducing length")
    ("ductSeg.absMaxTry",po::value<uint>()->default_value(10000),"max number of trial segments before completely giving up")
    // check if this needs to be changed
    ("ductSeg.roiStep",po::value<double>()->default_value(0.1),"step size for checking segment validity")
    ("ductSeg.segFrac",po::value<double>()->default_value(0.25),"fraction of branch length per segment")
    ;
  
  po::options_description vesselTreeOpt("Vessel tree options");
  vesselTreeOpt.add_options()
    ("vesselTree.maxBranch",po::value<uint>()->default_value(100),"Maximum number of branches")
    ("vesselTree.maxGen",po::value<uint>()->default_value(15),"Maximum generation")
    ("vesselTree.baseLength",po::value<double>()->default_value(12.0),"main branch length (mm)")
    ("vesselTree.initRad",po::value<double>()->default_value(2.0),"tree start radius")
    ("vesselTree.nFillX",po::value<uint>()->default_value(100),"number x voxels for density map")
    ("vesselTree.nFillY",po::value<uint>()->default_value(100),"number y voxels for density map")
    ("vesselTree.nFillZ",po::value<uint>()->default_value(100),"number z voxels for density map")
    ;

  po::options_description vesselBrOpt("Vessel branch options");
  vesselBrOpt.add_options()
    ("vesselBr.childMinRad",po::value<double>()->default_value(0.25),"min radius to have children")
    ("vesselBr.minRadFrac",po::value<double>()->default_value(0.6),"min start radius as fraction of parent radius")
    ("vesselBr.maxRadFrac",po::value<double>()->default_value(0.85),"max start radius as fraction of parent radius")
    ("vesselBr.lenShrink",po::value<double>()->default_value(0.5),"mean vessel length relative to parent")
    ("vesselBr.lenRange",po::value<double>()->default_value(0.1),"vessel length fraction variability")
    ("vesselBr.rotateJitter",po::value<double>()->default_value(0.1),"jitter in azimuthal angle of second child")

    ;
  
  po::options_description vesselSegOpt("Vessel segment options");
  vesselSegOpt.add_options()
    ("vesselSeg.radiusBetaA",po::value<double>()->default_value(2.0),"radius distribution shape parameter A")
    ("vesselSeg.radiusBetaB",po::value<double>()->default_value(2.0),"radius distribution shape parameter B")
    ("vesselSeg.maxCurvRad",po::value<double>()->default_value(20.0),"max radius of curvature")
    ("vesselSeg.maxCurvFrac",po::value<double>()->default_value(0.33),"max fraction of circle")
    ("vesselSeg.minEndRad",po::value<double>()->default_value(0.85),"min end radius as fraction of start radius")
    ("vesselSeg.maxEndRad",po::value<double>()->default_value(1.05),"max end radius as fraction of start radius")
    ("vesselSeg.angleWt",po::value<double>()->default_value(1.0),"cost function preferential angle weighting")
    ("vesselSeg.densityWt",po::value<double>()->default_value(5e-5),"cost function density weighting")
    ("vesselSeg.dirWt",po::value<double>()->default_value(5e-5),"cost function direction weighting")
    ("vesselSeg.numTry",po::value<uint>()->default_value(10),"number of trial segments")
    ("vesselSeg.maxTry",po::value<uint>()->default_value(100),"max number of trial segments before reducing length")
    ("vesselSeg.absMaxTry",po::value<uint>()->default_value(10000),"max number of trial segments before completely giving up")
    // check if this needs to be changed
    ("vesselSeg.roiStep",po::value<double>()->default_value(0.1),"step size for checking segment validity")
    ("vesselSeg.segFrac",po::value<double>()->default_value(0.25),"fraction of branch length per segment")
    ;

  // config file options
  po::options_description configFileOpt("Configuration file options");
  configFileOpt.add(baseOpt).add(shapeOpt);
  configFileOpt.add(ductTreeOpt).add(ductBrOpt).add(ductSegOpt);
  configFileOpt.add(vesselTreeOpt).add(vesselBrOpt).add(vesselSegOpt);
  configFileOpt.add(compartOpt).add(TDLUOpt).add(fatOpt);
  configFileOpt.add(voronOpt).add(perlinOpt).add(boundaryOpt);
  configFileOpt.add(bufferOpt).add(perturbOpt).add(ligOpt);

  // all of the options
  po::options_description all("All options");
  all.add_options()
    ("config,c", po::value<std::string>()->required(), "name of configuration file")
    ;
  all.add(configFileOpt);
  
  po::variables_map vm;
  
  // get configuration filename from command line
  po::store(parse_command_line(argc,argv,all), vm);
  std::string configFile = vm["config"].as<std::string>();
  
  // read configuration file
  ifstream inConfig(configFile.c_str());
  if(!inConfig){
    cout << "Can not open configuration file: " << configFile << "\n";
    return(1);
  } else {
    po::store(parse_config_file(inConfig, configFileOpt), vm);
    inConfig.close();
  };

  // load resolution variables

  double ures = vm["shape.ures"].as<double>();	// spacing in u-v space
  double vres = vm["shape.vres"].as<double>();
  double pointSep = vm["shape.pointSep"].as<double>(); // remove close points
  double featureAngle = vm["shape.featureAngle"].as<double>(); // angle to preserve when smoothing
  double targetReduction = vm["shape.targetReduction"].as<double>(); // fraction of triangles to decimate
  double imgRes = vm["base.imgRes"].as<double>();	// size of image voxels (mm)
  double skinThick = vm["base.skinThick"].as<double>(); // thickness of skin (mm)
  //std::string baseShapeFilename = "baseShape.vtp";
  std::string breastVoxelFilename = "breastVoxel.vti";
  //std::string breastPointsFilename = "breastPoints.vtp";
  //std::string segmentedVoxelFilename = "segmentedVoxel.vti";
  //std::string ductVoxelFilename = "final.vti";
  //std::string skinMeshFilename = "skinMesh.vtp";
  char outVTIFilename[128];
  //char outrawFilename[128];
  char outhdrFilename[128];
  char outgzFilename[128];


  double scaleFactor = 35.0;	// scale voxel size to millimeters (was 50)

  double nippleLen = vm["base.nippleLen"].as<double>();	// length of nipple (mm)
  double nippleRad = vm["base.nippleRad"].as<double>(); // radius of nipple (mm)
  double areolaRad = vm["base.areolaRad"].as<double>(); // radius of areola (mm)

  double voxelVol = imgRes*imgRes*imgRes;	// volume of a voxel (cubic mm)
	
  tissueStruct tissue;
	
  tissue.bg = 0;
  tissue.skin = 2;
  tissue.nipple = 33;
  tissue.fat = 1;
  tissue.cooper = 88;
  tissue.gland = 29;
  tissue.TDLU = 95;
  tissue.duct = 125;
  tissue.artery = 150;
  tissue.vein = 225 ;
  tissue.muscle = 40;

  const unsigned char glandBuffer = 222;
  const unsigned char innerVal = 175; // not needed after compartment creation
  const unsigned char boundVal = 50;  // edge of inner breast volume, not needed after skin generation
  
  const unsigned char compartmentVal[20] = {12, 13, 14, 15, 16,
					    17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
  
  const unsigned char compMax = 28;	// max compartment value for cooper's ligament part
  const unsigned char compMin = 12;

  // unligamented tissue classes
  const unsigned char ufat = 60;
  const unsigned char ugland = 61;
  const unsigned char uTDLU = 62;
  const unsigned char uduct = 63;

  // random number generator seed
  // use seed specified in config file, otherwise random seed
  int randSeed;
  if(vm.count("base.seed")){
    // seed specified in configuration file
    randSeed = vm["base.seed"].as<unsigned int>();
  } else {
    // generate random seed
    FILE *randFile;
    randFile = fopen("/dev/urandom","r");
    fread((char*)(&randSeed), sizeof(int),1,randFile);
    fclose(randFile);
  }

  // output base directory
  std::string outputDir = vm["base.outputDir"].as<std::string>();
  
  // exists?
  if(access(outputDir.c_str(),F_OK)){
    // directory doesn't exist
    // try to create it
    if(mkdir(outputDir.c_str(),S_IRWXU) == -1){
      // couldn't create directory
      cerr << "Could not create directory " << outputDir << "\n";
      cerr << "Exiting...\n";
      return(1);
    }
  }

  // is it a directory?
  struct stat fileInfo;
  
  stat(outputDir.c_str(), &fileInfo);

  if(S_ISDIR(fileInfo.st_mode)){
    // it's a directory
    // make it the working directory
    if(chdir(outputDir.c_str()) == -1){
      cerr << "Could not access directory " << outputDir << "\n";
      cerr << "Exiting...\n";
      return(1);
    }
  } else {
    // error, not a directory
    cerr << "Specified path is not a valid directory\n";
    cerr << "Exiting...\n";
    return(1);
  }

  // copy over config file
  
  char cfgOutFilename[128];
  sprintf(cfgOutFilename,"%s/p_%d.cfg", outputDir.c_str(),randSeed);
	
  FILE *cfgCopy = fopen(cfgOutFilename,"wt");
  FILE *cfgRead = fopen(configFile.c_str(),"rt");
  
  int readChar = getc(cfgRead);
  while(readChar != EOF){
    putc(readChar, cfgCopy);
    readChar = getc(cfgRead);
  }
  fclose(cfgCopy);
  fclose(cfgRead);

  sprintf(outVTIFilename,"%s/p_%d.vti", outputDir.c_str(),randSeed);
  //sprintf(outrawFilename,"%s/p_%d.zraw", outputDir.c_str(),randSeed);
  sprintf(outhdrFilename,"%s/p_%d.mhd", outputDir.c_str(),randSeed);
  sprintf(outgzFilename,"%s/p_%d.raw.gz", outputDir.c_str(),randSeed);

  // shape parameters

  // base shape coefficients
  double a1b =  vm["shape.a1b"].as<double>();
  double a1t = vm["shape.a1t"].as<double>();
  double a2l = vm["shape.a2l"].as<double>();
  double a2r = vm["shape.a2r"].as<double>();
  double a3 = vm["shape.a3"].as<double>();
  double eps1 = vm["shape.eps1"].as<double>();
  double eps2 = vm["shape.eps2"].as<double>();

  // breast side
  bool leftSide = vm["base.leftBreast"].as<bool>();

  // ptosis parameters
  bool doPtosis = vm["shape.doPtosis"].as<bool>();
  double b0 = vm["shape.ptosisB0"].as<double>();
  double b1 = vm["shape.ptosisB1"].as<double>();

  // turn parameters
  bool doTurn = vm["shape.doTurn"].as<bool>();
  double c0 = vm["shape.turnC0"].as<double>();
  double c1 = vm["shape.turnC1"].as<double>();

  // top shape parameters
  bool doTopShape = vm["shape.doTopShape"].as<bool>();
  double s0 = vm["shape.topShapeS0"].as<double>();
  double t0 = vm["shape.topShapeT0"].as<double>();
  double s1 = vm["shape.topShapeS1"].as<double>();
  double t1 = vm["shape.topShapeT1"].as<double>();

  // derived top shape parameters
  double At = -0.5*t0-3.0*s0-3.0*s1+0.5*t1;
  double Bt = 1.5*t0+8.0*s0+7.0*s1-t1;
  double Ct = -1.5*t0-6.0*s0-4.0*s1+0.5*t1;
  double Dt = 0.5*t0;
  double Et = s0;
  double Ft = 1.0;

  // flatten side parameters
  bool doFlattenSide = vm["shape.doFlattenSide"].as<bool>();
  double g0 = vm["shape.flattenSideG0"].as<double>();
  double g1 = vm["shape.flattenSideG1"].as<double>();

  // derived parameters for flatten side
  double Af = g1+2.0-2.0*g0;
  double Bf = -g1-3.0+3.0*g0;
  double Cf = 0.0;
  double Df = 1.0;

  // turn top parameters
  bool doTurnTop = vm["shape.doTurnTop"].as<bool>();
  double h0 = vm["shape.turnTopH0"].as<double>();
  double h1 = vm["shape.turnTopH1"].as<double>();

  // start a random number generator
  vtkSmartPointer<vtkMinimalStandardRandomSequence> rgen =
    vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();

  rgen->SetSeed(randSeed);

  /***********************
	Shape
  ***********************/

  // create base shape
  //cout << "Initializing point cloud...";
  // point positions
  double uval,vval,xval,yval,zval;
  vtkIdType myId;
  
  vtkSmartPointer<vtkPoints> brightFront =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> brightBack =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> brightRing =
    vtkSmartPointer<vtkPoints>::New();
		
  vtkSmartPointer<vtkPoints> trightFront =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> trightBack =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> trightRing =
    vtkSmartPointer<vtkPoints>::New();

  vtkSmartPointer<vtkPoints> tleftFront =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> tleftBack =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> tleftRing =
    vtkSmartPointer<vtkPoints>::New();
	
  vtkSmartPointer<vtkPoints> bleftFront =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> bleftBack =
    vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> bleftRing =
    vtkSmartPointer<vtkPoints>::New();
	
  vtkIdType centerPt;
	
#pragma omp sections private(uval,vval,xval,yval,zval,myId)
  {
#pragma omp section
    {
      // bottom right breast
      vval = -1.0*pi + vres;
      while(vval <= -0.5*pi){
	uval = 0.0;
	// ring is uval == 0
	xval = 0.0;
	yval = pow(a2r,eps1)*pow(sin(vval),eps2);
	zval = pow(a1b,eps1)*pow(cos(vval),eps2);
	myId = brightRing->InsertNextPoint(xval,yval,zval);
	uval += ures;
			
	while(uval <= 0.5*pi){
	  xval = pow(a3*sin(uval),eps1);
	  yval = pow(a2r*cos(uval),eps1)*pow(sin(vval),eps2);
	  zval = pow(a1b*cos(uval),eps1)*pow(cos(vval),eps2);
	  // add point to front
	  myId = brightFront->InsertNextPoint(xval,yval,zval);
	  // add backplane point
	  myId = brightBack->InsertNextPoint(0.0,yval,zval);
	  uval += ures;
	}
	vval += vres;
      }
	
      // add center point to bottom-right breast front and back and skip
      centerPt = brightFront->InsertNextPoint(pow(a3,eps1),0.0,0.0);
      centerPt = brightBack->InsertNextPoint(0.0,0.0,0.0);
    }
    
#pragma omp section
    {
      // top right breast
      
      vval = -0.5*pi + vres;
      while(vval <= 0.0){
	uval = 0.0;
	xval = 0.0;
	yval = pow(a2r,eps1)*pow(sin(vval),eps2);
	zval = pow(a1t,eps1)*pow(cos(vval),eps2);
	myId = trightRing->InsertNextPoint(xval,yval,zval);
	uval += ures;
			
	while(uval <= 0.5*pi){
	  xval = pow(a3*sin(uval),eps1);
	  yval = pow(a2r*cos(uval),eps1)*pow(sin(vval),eps2);
	  zval = pow(a1t*cos(uval),eps1)*pow(cos(vval),eps2);
	  // add the point
	  myId = trightFront->InsertNextPoint(xval,yval,zval);
	  // add backplane point
	  myId = trightBack->InsertNextPoint(0.0,yval,zval);
	  uval += ures;
	}
	vval += vres;
      }
    }
		
#pragma omp section
    {
      // top left breast
		
      vval = 0.0 + vres;
      while(vval <= 0.5*pi){
	uval = 0.0;
	xval = 0.0;
	yval = pow(a2l,eps1)*pow(sin(vval),eps2);
	zval = pow(a1t,eps1)*pow(cos(vval),eps2);
	myId = tleftRing->InsertNextPoint(xval,yval,zval);
	uval += ures;
	
	while(uval <= 0.5*pi){
	  xval = pow(a3*sin(uval),eps1);
	  yval = pow(a2l*cos(uval),eps1)*pow(sin(vval),eps2);
	  zval = pow(a1t*cos(uval),eps1)*pow(cos(vval),eps2);
	  // add the point
	  myId = tleftFront->InsertNextPoint(xval,yval,zval);
	  // add backplane point
	  myId = tleftBack->InsertNextPoint(0.0,yval,zval);
	  uval += ures;
	}
	vval += vres;
      }
    }

#pragma omp section
    {
      // bottom left breast
      
      vval = 0.5*pi + vres;
      while(vval <= pi){
	uval = 0.0;
	xval = 0.0;
	yval = pow(a2l,eps1)*pow(sin(vval),eps2);
	zval = pow(a1b,eps1)*pow(cos(vval),eps2);
	myId = bleftRing->InsertNextPoint(xval,yval,zval);
	uval += ures;
			
	while(uval <= 0.5*pi){
	  xval = pow(a3*sin(uval),eps1);
	  yval = pow(a2l*cos(uval),eps1)*pow(sin(vval),eps2);
	  zval = pow(a1b*cos(uval),eps1)*pow(cos(vval),eps2);
	  // add the point
	  myId = bleftFront->InsertNextPoint(xval,yval,zval);
	  // add backplane point
	  myId = bleftBack->InsertNextPoint(0.0,yval,zval);
	  uval += ures;
	}
	vval += vres;
      }
    }

  }  // end omp sections
	
  //cout << "done.\n";

  // have base shape
  // do deformations

  //cout << "Doing deformation...";

  // top shape
  if(doTopShape){
    vtkIdType npts = trightFront->GetNumberOfPoints();
    uval = ures; // need original u value
    for(int i=0; i<npts; i++){
      double pval[3];
      double u2 = uval*2.0/pi;
      // front
      trightFront->GetPoint(i, pval);
      pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
			 Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
      trightFront->SetPoint(i,pval);
      // back
      trightBack->GetPoint(i, pval);
      pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
			 Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
      trightBack->SetPoint(i,pval);
      
      uval += ures;
      if(uval > 0.5*pi){
	uval = ures;
      }
    }
		
    // ring
    npts = trightRing->GetNumberOfPoints();
    for(int i=0; i<npts; i++){
      double pval[3];
      trightRing->GetPoint(i, pval);
      pval[2] = pval[2]*Ft;
      trightRing->SetPoint(i,pval);
    }
		
    npts = tleftFront->GetNumberOfPoints();
    uval = ures; // need original u value
    for(int i=0; i<npts; i++){
      double pval[3];
      double u2 = uval*2.0/pi;
      // front
      tleftFront->GetPoint(i, pval);
      pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
			 Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
      tleftFront->SetPoint(i,pval);
      // back
      tleftBack->GetPoint(i, pval);
      pval[2] = pval[2]*(At*pow(u2,5.0)+Bt*pow(u2,4.0)+
			 Ct*pow(u2,3.0)+Dt*u2*u2+Et*u2+Ft);
      tleftBack->SetPoint(i,pval);
			
      uval += ures;
      if(uval > 0.5*pi){
	uval = ures;
      }
    }
		
    // ring
    npts = tleftRing->GetNumberOfPoints();
    for(int i=0; i<npts; i++){
      double pval[3];
      tleftRing->GetPoint(i, pval);
      pval[2] = pval[2]*Ft;
      tleftRing->SetPoint(i,pval);
    }				
  }
	
  // flatten side
  if(doFlattenSide){
    if(leftSide){
      trightFront->ComputeBounds();
      double bound[6];
      trightFront->GetBounds(bound);
      double scale;
      if(fabs(bound[2]) > fabs(bound[3])){
	scale = fabs(bound[2]);
      } else {
	scale = fabs(bound[3]);
      }
      vtkIdType npts = trightFront->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	// front
	trightFront->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	trightFront->SetPoint(i,pval);
	// back
	trightBack->GetPoint(i,pval);
	yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	trightBack->SetPoint(i,pval);
      }
      // ring
      npts = trightRing->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	trightRing->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	trightRing->SetPoint(i,pval);
      }
						
      brightFront->ComputeBounds();
      brightFront->GetBounds(bound);
      if(fabs(bound[2]) > fabs(bound[3])){
	scale = fabs(bound[2]);
      } else {
	scale = fabs(bound[3]);
      }
      npts = brightFront->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	//front
	brightFront->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	brightFront->SetPoint(i,pval);
	// back
	brightBack->GetPoint(i,pval);
	yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	brightBack->SetPoint(i,pval);
      }
      // ring
      npts = brightRing->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	brightRing->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	brightRing->SetPoint(i,pval);
      }
    } else {
      // right side
      tleftFront->ComputeBounds();
      double bound[6];
      tleftFront->GetBounds(bound);
      double scale;
      if(fabs(bound[2]) > fabs(bound[3])){
	scale = fabs(bound[2]);
      } else {
	scale = fabs(bound[3]);
      }
      vtkIdType npts = tleftFront->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	// front
	tleftFront->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	tleftFront->SetPoint(i,pval);
	// back
	tleftBack->GetPoint(i,pval);
	yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	tleftBack->SetPoint(i,pval);
      }
      // ring
      npts = tleftRing->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	tleftRing->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	tleftRing->SetPoint(i,pval);
      }
      
      bleftFront->ComputeBounds();
      bleftFront->GetBounds(bound);
      if(fabs(bound[2]) > fabs(bound[3])){
	scale = fabs(bound[2]);
      } else {
	scale = fabs(bound[3]);
      }
      npts = bleftFront->GetNumberOfPoints();
#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	// front
	bleftFront->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	bleftFront->SetPoint(i,pval);
	// back
	bleftBack->GetPoint(i,pval);
	yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	bleftBack->SetPoint(i,pval);
      }
      // ring
      npts = bleftRing->GetNumberOfPoints();
      //#pragma omp parallel for
      for(int i=0; i<npts; i++){
	double pval[3];
	bleftRing->GetPoint(i,pval);
	double yv = fabs(pval[1]/scale);
	pval[1] = pval[1]*(Af*yv*yv*yv + Bf*yv*yv +
			   Cf*yv + Df);
	bleftRing->SetPoint(i,pval);
      }
    }
  }
	
  // turn top
  if(doTurnTop){
    trightFront->ComputeBounds();
    double bound[6];
    trightFront->GetBounds(bound);
    double scale;
    if(fabs(bound[4]) > fabs(bound[5])){
      scale = fabs(bound[4]);
    } else {
      scale = fabs(bound[5]);
    }
    vtkIdType npts = trightFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      trightFront->GetPoint(i,pval);
      pval[1] = pval[1] - h0*pval[2]/scale -
	h1*pval[2]*pval[2]/scale/scale;
      trightFront->SetPoint(i,pval);
      // back
      trightBack->GetPoint(i,pval);
      pval[1] = pval[1] - h0*pval[2]/scale -
	h1*pval[2]*pval[2]/scale/scale;
      trightBack->SetPoint(i,pval);
    }
    // ring
    npts = trightRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      trightRing->GetPoint(i,pval);
      pval[1] = pval[1] - h0*pval[2]/scale -
	h1*pval[2]*pval[2]/scale/scale;
      trightRing->SetPoint(i,pval);
    }
		
    tleftFront->ComputeBounds();
    tleftFront->GetBounds(bound);
    if(fabs(bound[4]) > fabs(bound[5])){
      scale = fabs(bound[4]);
    } else {
      scale = fabs(bound[5]);
    }
    npts = tleftFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      tleftFront->GetPoint(i,pval);
      pval[1] = pval[1] - h0*pval[2]/scale -
	h1*pval[2]*pval[2]/scale/scale;
      tleftFront->SetPoint(i,pval);
      // back
      tleftBack->GetPoint(i,pval);
      pval[1] = pval[1] - h0*pval[2]/scale -
	h1*pval[2]*pval[2]/scale/scale;
      tleftBack->SetPoint(i,pval);
    }
    // ring
    npts = tleftRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      tleftRing->GetPoint(i,pval);
      pval[1] = pval[1] - h0*pval[2]/scale -
	h1*pval[2]*pval[2]/scale/scale;
      tleftRing->SetPoint(i,pval);
    }
  }

  // ptosis
  if(doPtosis){
    vtkIdType npts = trightFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      trightFront->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      trightFront->SetPoint(i,pval);
      // back
      trightBack->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      trightBack->SetPoint(i,pval);
    }
    // ring
    npts = trightRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      trightRing->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      trightRing->SetPoint(i,pval);
    }
    
    npts = tleftFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      tleftFront->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      tleftFront->SetPoint(i,pval);
      // back
      tleftBack->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      tleftBack->SetPoint(i,pval);
    }
    // ring
    npts = tleftRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      tleftRing->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      tleftRing->SetPoint(i,pval);
    }
    
    npts = brightFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      brightFront->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      brightFront->SetPoint(i,pval);
      // back
      brightBack->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      brightBack->SetPoint(i,pval);
    }
    // ring
    npts = brightRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      brightRing->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      brightRing->SetPoint(i,pval);
    }
		
    npts = bleftFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      bleftFront->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      bleftFront->SetPoint(i,pval);
      // back
      bleftBack->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      bleftBack->SetPoint(i,pval);
    }
    // ring
    npts = bleftRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      bleftRing->GetPoint(i,pval);
      pval[2] = pval[2] - (b0*pval[0] + b1*pval[0]*pval[0]);
      bleftRing->SetPoint(i,pval);
    }
  }

  // turn
  if(doTurn){
    vtkIdType npts = trightFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      trightFront->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      trightFront->SetPoint(i,pval);
      // back
      trightBack->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      trightBack->SetPoint(i,pval);
    }
    // ring
    npts = trightRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      trightRing->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      trightRing->SetPoint(i,pval);
    }
    
    npts = tleftFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      tleftFront->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      tleftFront->SetPoint(i,pval);
      // back
      tleftBack->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      tleftBack->SetPoint(i,pval);
    }
    // ring
    npts = tleftRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      tleftRing->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      tleftRing->SetPoint(i,pval);
    }
		
    npts = brightFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      brightFront->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      brightFront->SetPoint(i,pval);
      // back
      brightBack->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      brightBack->SetPoint(i,pval);
    }
    // ring
    npts = brightRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      brightRing->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      brightRing->SetPoint(i,pval);
    }
    
    npts = bleftFront->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      // front
      bleftFront->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      bleftFront->SetPoint(i,pval);
      // back
      bleftBack->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      bleftBack->SetPoint(i,pval);
    }
    // ring
    npts = bleftRing->GetNumberOfPoints();
#pragma omp parallel for
    for(int i=0; i<npts; i++){
      double pval[3];
      bleftRing->GetPoint(i,pval);
      pval[1] = pval[1] + (c0*pval[0] + c1*pval[0]*pval[0]);
      bleftRing->SetPoint(i,pval);
    }
  }

  // get nipple position
  double nipplePos[3];
  brightFront->GetPoint(centerPt,nipplePos);

  // scale nipple position
  for(int i=0; i<3; i++){
    nipplePos[i] = nipplePos[i]*scaleFactor;
  }

  // aggregate front points
  vtkSmartPointer<vtkPoints> frontPts =
    vtkSmartPointer<vtkPoints>::New();
	
  for(vtkIdType i=0; i<brightFront->GetNumberOfPoints(); i++){
    double t[3];
    brightFront->GetPoint(i, t);
    frontPts->InsertNextPoint(t);
  }
  for(vtkIdType i=0; i<trightFront->GetNumberOfPoints(); i++){
    double t[3];
    trightFront->GetPoint(i, t);
    frontPts->InsertNextPoint(t);
  }
  for(vtkIdType i=0; i<bleftFront->GetNumberOfPoints(); i++){
    double t[3];
    bleftFront->GetPoint(i, t);
    frontPts->InsertNextPoint(t);
  }
  for(vtkIdType i=0; i<tleftFront->GetNumberOfPoints(); i++){
    double t[3];
    tleftFront->GetPoint(i, t);
    frontPts->InsertNextPoint(t);
  }
  vtkIdType numFrontPts = frontPts->GetNumberOfPoints();
	
  // decimate front points
  vtkSmartPointer<vtkPolyData> frontPoly = 
    vtkSmartPointer<vtkPolyData>::New();
  frontPoly->SetPoints(frontPts);
  vtkSmartPointer<vtkVertexGlyphFilter> frontVertAdd =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  frontVertAdd->SetInput(frontPoly);
#else
  frontVertAdd->SetInputData(frontPoly);
#endif

  frontVertAdd->Update();
  vtkSmartPointer<vtkCleanPolyData> cleanFront =
    vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
  cleanFront->SetInput(frontVertAdd->GetOutput());
#else
  cleanFront->SetInputConnection(frontVertAdd->GetOutputPort());
#endif
  cleanFront->SetTolerance(pointSep);
  cleanFront->Update();	
	
  // aggregate and shift back points
  vtkSmartPointer<vtkPoints> backPts =
    vtkSmartPointer<vtkPoints>::New();
	
  // ringWidth and ringSep converted to non-physical coordinates
  double ringWidth = vm["shape.ringWidth"].as<double>();
  double ringSep = vm["shape.ringSep"].as<double>();
  double ringWidthOrig = ringWidth/scaleFactor;
  double ringSepOrig = ringSep/scaleFactor;
  double backPos = (floor(ringWidthOrig/ringSepOrig)+1)*ringSepOrig;
  
  for(vtkIdType i=0; i<brightBack->GetNumberOfPoints(); i++){
    double t[3];
    brightBack->GetPoint(i, t);
    t[0] = -backPos;
    backPts->InsertNextPoint(t);
  }
  for(vtkIdType i=0; i<trightBack->GetNumberOfPoints(); i++){
    double t[3];
    trightBack->GetPoint(i, t);
    t[0] = -backPos;
    backPts->InsertNextPoint(t);
  }
  for(vtkIdType i=0; i<bleftBack->GetNumberOfPoints(); i++){
    double t[3];
    bleftBack->GetPoint(i, t);
    t[0] = -backPos;
    backPts->InsertNextPoint(t);
  }
  for(vtkIdType i=0; i<tleftBack->GetNumberOfPoints(); i++){
    double t[3];
    tleftBack->GetPoint(i, t);
    t[0] = -backPos;
    backPts->InsertNextPoint(t);
  }
  vtkIdType numBackPts = backPts->GetNumberOfPoints();
  
  // decimate back points
  vtkSmartPointer<vtkPolyData> backPoly = 
    vtkSmartPointer<vtkPolyData>::New();
  backPoly->SetPoints(backPts);
  vtkSmartPointer<vtkVertexGlyphFilter> backVertAdd =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  backVertAdd->SetInput(backPoly);
#else
  backVertAdd->SetInputData(backPoly);
#endif
  backVertAdd->Update();
  vtkSmartPointer<vtkCleanPolyData> cleanBack =
    vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
  cleanBack->SetInput(backVertAdd->GetOutput());
#else
  cleanBack->SetInputConnection(backVertAdd->GetOutputPort());
#endif
  cleanBack->SetTolerance(pointSep);
  cleanBack->Update();	
	
  // aggregate ring points and create thickness
  vtkSmartPointer<vtkPoints> ringPts =
    vtkSmartPointer<vtkPoints>::New();
	
  for(vtkIdType i=0; i<brightRing->GetNumberOfPoints(); i++){
    double t[3];
    brightRing->GetPoint(i, t);
    ringPts->InsertNextPoint(t);
    double xstep = ringSepOrig;
    while(xstep < ringWidthOrig){
      t[0] = -xstep;
      ringPts->InsertNextPoint(t);
      xstep += ringSepOrig;
    }
  }
  for(vtkIdType i=0; i<trightRing->GetNumberOfPoints(); i++){
    double t[3];
    trightRing->GetPoint(i, t);
    ringPts->InsertNextPoint(t);
    double xstep = ringSepOrig;
    while(xstep < ringWidthOrig){
      t[0] = -xstep;
      ringPts->InsertNextPoint(t);
      xstep += ringSepOrig;
    }
  }
  for(vtkIdType i=0; i<bleftRing->GetNumberOfPoints(); i++){
    double t[3];
    bleftRing->GetPoint(i, t);
    ringPts->InsertNextPoint(t);
    double xstep = ringSepOrig;
    while(xstep < ringWidthOrig){
      t[0] = -xstep;
      ringPts->InsertNextPoint(t);
      xstep += ringSepOrig;
    }
  }
  for(vtkIdType i=0; i<tleftRing->GetNumberOfPoints(); i++){
    double t[3];
    tleftRing->GetPoint(i, t);
    ringPts->InsertNextPoint(t);
    double xstep = ringSepOrig;
    while(xstep < ringWidthOrig){
      t[0] = -xstep;
      ringPts->InsertNextPoint(t);
      xstep += ringSepOrig;
    }
  }
  vtkIdType numRingPts = ringPts->GetNumberOfPoints();
  
  // decimate ring points
  vtkSmartPointer<vtkPolyData> ringPoly = 
    vtkSmartPointer<vtkPolyData>::New();
  ringPoly->SetPoints(ringPts);
  vtkSmartPointer<vtkVertexGlyphFilter> ringVertAdd =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  ringVertAdd->SetInput(ringPoly);
#else
  ringVertAdd->SetInputData(ringPoly);
#endif
  ringVertAdd->Update();
  vtkSmartPointer<vtkCleanPolyData> cleanRing =
    vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
  cleanRing->SetInput(ringVertAdd->GetOutput());
#else
  cleanRing->SetInputConnection(ringVertAdd->GetOutputPort());
#endif
  cleanRing->SetTolerance(pointSep);
  cleanRing->Update();	
	
  // allocate space for breast points
  vtkSmartPointer<vtkPoints> breastPts =
    vtkSmartPointer<vtkPoints>::New();

  vtkIdType numBreastPts = cleanFront->GetOutput()->GetNumberOfPoints();
  numBreastPts += cleanBack->GetOutput()->GetNumberOfPoints();
  numBreastPts += cleanRing->GetOutput()->GetNumberOfPoints();
  breastPts->SetNumberOfPoints(numBreastPts);
  
  // cell array for verticies
  vtkSmartPointer<vtkCellArray> breastVerts =
    vtkSmartPointer<vtkCellArray>::New();
	
  // scale to physical units
  vtkIdType totalCount = 0;
  vtkIdType currentPoints = cleanFront->GetOutput()->GetNumberOfPoints();

  for(vtkIdType i=0; i<currentPoints; i++){
    double t[3];
		
    cleanFront->GetOutput()->GetPoint(i,t);
    for(int j=0; j<3; j++){
      t[j] = t[j]*scaleFactor;
    }
    breastPts->SetPoint(totalCount,t);
    breastVerts->InsertNextCell(1,&totalCount); 
    totalCount++;
  }

  // calculate scaled back points
  currentPoints = cleanBack->GetOutput()->GetNumberOfPoints();
  for(vtkIdType i=0; i<currentPoints; i++){
    double t[3];
    cleanBack->GetOutput()->GetPoint(i,t);
    for(int j=0; j<3; j++){
      t[j] = t[j]*scaleFactor;
    }
    breastPts->SetPoint(totalCount,t);
    breastVerts->InsertNextCell(1,&totalCount); 
    totalCount++;
  }
	
  // calculate scaled ring points and normals
  currentPoints = cleanRing->GetOutput()->GetNumberOfPoints();

  for(vtkIdType i=0; i<currentPoints; i++){
    double t[3];
    cleanRing->GetOutput()->GetPoint(i,t);
    for(int j=0; j<3; j++){
      t[j] = t[j]*scaleFactor;
    }
    breastPts->SetPoint(totalCount,t);
    breastVerts->InsertNextCell(1,&totalCount); 
    totalCount++;
  }

  //cout << "done.\n";

  //cout << "Meshing base shape...";

  // create polydata
  vtkSmartPointer<vtkPolyData> breastPoly =
    vtkSmartPointer<vtkPolyData>::New();
  breastPoly->SetPoints(breastPts);
	
  // add verticies to polydata
  breastPoly->SetVerts(breastVerts);
	
  // debug
  //vtkSmartPointer<vtkXMLPolyDataWriter> writerd =
  //	vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //writerd->SetFileName(breastPointsFilename.c_str());
  //#if VTK_MAJOR_VERSION <= 5
  //writerd->SetInput(breastPoly);
  //#else
  //writerd->SetInputData(breastPoly);
  //#endif
  //writerd->Write();

  // surface mesh reconstruction
  vtkSmartPointer<vtkSurfaceReconstructionFilter> breastFilter =
    vtkSmartPointer<vtkSurfaceReconstructionFilter>::New();
		
#if VTK_MAJOR_VERSION <= 5
  breastFilter->SetInput(breastPoly);
#else
  breastFilter->SetInputData(breastPoly);
#endif
  breastFilter->Update();
	
  // fix bug in vtkSurfaceReconstructionFilter
  double pointSetBounds[6];
  breastPts->GetBounds(pointSetBounds);
	
  vtkSmartPointer<vtkContourFilter> breastContourFilter =
    vtkSmartPointer<vtkContourFilter>::New();
		
#if VTK_MAJOR_VERSION <= 5
  breastContourFilter->SetInput(breastFilter->GetOutput());
#else
  breastContourFilter->SetInputConnection(breastFilter->GetOutputPort());
#endif
  breastContourFilter->ComputeNormalsOn();
  breastContourFilter->SetValue(0, 0.0);
  breastContourFilter->Update();
	
  vtkSmartPointer<vtkReverseSense> breastSenseFilter =
    vtkSmartPointer<vtkReverseSense>::New();
	
#if VTK_MAJOR_VERSION <= 5
  breastSenseFilter->SetInput(breastContourFilter->GetOutput());
#else
  breastSenseFilter->SetInputConnection(breastContourFilter->GetOutputPort());
#endif
  breastSenseFilter->ReverseCellsOn();
  breastSenseFilter->ReverseNormalsOn();
  breastSenseFilter->Update();

  double surfaceBounds[6];
  breastSenseFilter->GetOutput()->GetBounds(surfaceBounds);
	
  double scaling[3];
	
  scaling[0] = (pointSetBounds[1]-pointSetBounds[0])/(surfaceBounds[1]-surfaceBounds[0]);
  scaling[1] = (pointSetBounds[3]-pointSetBounds[2])/(surfaceBounds[3]-surfaceBounds[2]);
  scaling[2] = (pointSetBounds[5]-pointSetBounds[4])/(surfaceBounds[5]-surfaceBounds[4]);

  vtkSmartPointer<vtkTransform> transp = 
    vtkSmartPointer<vtkTransform>::New();
  transp->Translate(pointSetBounds[0], pointSetBounds[2], pointSetBounds[4]);
  transp->Scale(scaling[0], scaling[1], scaling[2]);
  transp->Translate(-surfaceBounds[0], -surfaceBounds[2], -surfaceBounds[4]);

  vtkSmartPointer<vtkTransformPolyDataFilter> fixed = 
    vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  fixed->SetInputData(breastSenseFilter->GetOutput());
#if VTK_MAJOR_VERSION <= 5
  fixed->SetInput(breastSenseFilter->GetOutput());
#else
  fixed->SetInputData(breastSenseFilter->GetOutput());
#endif
  fixed->SetTransform(transp);
  fixed->Update();

  // clean mesh
  vtkSmartPointer<vtkCleanPolyData> cleanPoly =
    vtkSmartPointer<vtkCleanPolyData>::New();
#if VTK_MAJOR_VERSION <= 5
  cleanPoly->SetInput(fixed->GetOutput());
#else
  cleanPoly->SetInputConnection(fixed->GetOutputPort());
#endif
	
  // decimate mesh
  vtkSmartPointer<vtkDecimatePro> decimate = 
    vtkSmartPointer<vtkDecimatePro>::New();
  // debug
  decimate->SetTargetReduction(0.25);
#if VTK_MAJOR_VERSION <= 5
  decimate->SetInput(cleanPoly->GetOutput());
#else
  decimate->SetInputConnection(cleanPoly->GetOutputPort());
#endif
  decimate->Update();

  vtkSmartPointer<vtkPolyData> innerPoly =
    vtkSmartPointer<vtkPolyData>::New();
  innerPoly->ShallowCopy(decimate->GetOutput());

  // at this point have base surface
  // save base surface to disk
  //vtkSmartPointer<vtkXMLPolyDataWriter> writer1 =
  //	vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //writer1->SetFileName(baseShapeFilename.c_str());
//#if VTK_MAJOR_VERSION <= 5
  //writer1->SetInput(decimate->GetOutput());
//#else
  //writer1->SetInputConnection(decimate->GetOutputPort());
//#endif
  //writer1->Write();

  //cout << "done.\n";

  // timing
  //time(&currentTime);
  //cout << "Time interval 1: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  //cout << "Creating voxelization...";
  // convert inner breast to voxel representation

  // create 3d imagedata
  vtkSmartPointer<vtkImageData> breast =
    vtkSmartPointer<vtkImageData>::New();

  double baseBound[6];
  innerPoly->GetBounds(baseBound);
  
  // extend bounds slightly to make room for nipple
  double baseFOV[3];
  for(int i=0; i<3; i++){
    baseFOV[i] = fabs(baseBound[2*i+1] - baseBound[2*i]);
  }

  double finalBound[6];
  for(int i=0; i<6; i++){
    finalBound[i] = baseBound[i];
  }
  finalBound[0] = finalBound[0] + 1.0;	// move backplane in by 1 mm to avoid surface edge
  finalBound[1] = finalBound[1] + nippleLen + 2.0;	// make room for nipple
  finalBound[2] = finalBound[2] - 0.01*baseFOV[1];
  finalBound[3] = finalBound[3] + 0.01*baseFOV[1];
  finalBound[4] = finalBound[4] - 0.01*baseFOV[2];
  finalBound[5] = finalBound[5] + 0.01*baseFOV[2];

  double spacing[3];
  for(int i=0; i<3; i++){
    spacing[i] = imgRes;
  }
  breast->SetSpacing(spacing);

  int dim[3];
  for(int i=0; i<3; i++){
    dim[i] = static_cast<int>(ceil((finalBound[2*i+1]-finalBound[2*i])/spacing[i]));
  }

  breast->SetExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1);
  
  double origin[3];
  for(int i=0; i<3; i++){
    origin[i] = finalBound[2*i]+spacing[i]/2.0;
  }
  breast->SetOrigin(origin);

  // allocate unsigned char
#if VTK_MAJOR_VERSION <= 5
  breast->SetNumberOfScalarComponents(1);
  breast->SetScalarTypeToUnsignedChar();
  breast->AllocateScalars();
#else
  breast->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif

  int originIndex[3] = {0, 0, 0};
  double originCoords[3];
  breast->GetPoint(breast->ComputePointId(originIndex),originCoords);

  // initialize to zero
  unsigned char* voxVal = static_cast<unsigned char *>(breast->GetScalarPointer());
  const long long int numElements = dim[0]*dim[1]*dim[2];
  for(long long int i=0; i<numElements; i++){
    *voxVal = tissue.bg;
    voxVal++;
  }

  // voxelize

  // list of boundary voxels
  vtkSmartPointer<vtkIdList> boundaryList =
    vtkSmartPointer<vtkIdList>::New();

  vtkSmartPointer<vtkIdList> boundaryList1 =
    vtkSmartPointer<vtkIdList>::New();

  vtkSmartPointer<vtkIdList> boundaryList2 =
    vtkSmartPointer<vtkIdList>::New();

  vtkSmartPointer<vtkIdList> boundaryList3 =
    vtkSmartPointer<vtkIdList>::New();



  int maxThread = omp_get_max_threads();

  // boundary sub-lists
  vtkSmartPointer<vtkIdList> *subList1 = new vtkSmartPointer<vtkIdList>[maxThread];
  vtkSmartPointer<vtkIdList> *subList2 = new vtkSmartPointer<vtkIdList>[maxThread];
  vtkSmartPointer<vtkIdList> *subList3 = new vtkSmartPointer<vtkIdList>[maxThread];
  for(int i=0; i<maxThread; i++){
    subList1[i] = vtkSmartPointer<vtkIdList>::New();
    subList2[i] = vtkSmartPointer<vtkIdList>::New();
    subList3[i] = vtkSmartPointer<vtkIdList>::New();
  }

#pragma omp parallel num_threads(maxThread)
  { 
    int numThread = omp_get_num_threads();
    int myThread = omp_get_thread_num();

    vtkSmartPointer<vtkPolyData> myPoly =
      vtkSmartPointer<vtkPolyData>::New();
    myPoly->DeepCopy(innerPoly);
 
    // Create the tree
    vtkSmartPointer<vtkCellLocator> innerLocator =
      vtkSmartPointer<vtkCellLocator>::New();
    innerLocator->SetDataSet(myPoly);
    innerLocator->BuildLocator();

    // cout << "Thread " << myThread << ": locator built\n";

    // find intersect with top and bottom surface to voxelize breast
    // iterate over x and y values
    for(int i=myThread; i<dim[0]; i+=numThread){
      double xpos = origin[0]+i*spacing[0];
      int ijk[3];
      ijk[0] = i;
      for(int j=0; j<dim[1]; j++){
      
	double ypos = origin[1]+j*spacing[1];
	ijk[1] = j;
	// calculate z position of top surface
	double lineStart[3]; // end points of line
	double lineEnd[3];
      
	double tol = 0.005;
	double tval;
	vtkIdType intersectCell;
	int subId;

	double intersect[3]; // output position
	double pcoords[3];
      
	lineStart[0] = xpos;
	lineStart[1] = ypos;
	lineStart[2] = baseBound[4];

	lineEnd[0] = xpos;
	lineEnd[1] = ypos;
	lineEnd[2] = baseBound[5];

	if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					   tval, intersect, pcoords, subId, intersectCell)){

	  // found intersection
	  double topZ = intersect[2];

	  // do other direction
	  lineStart[2] = baseBound[5];
	  lineEnd[2] = baseBound[4];

	  if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					     tval, intersect, pcoords, subId, intersectCell)){

	    double bottomZ = intersect[2];
	    
	    // find nearest voxels to intersections
	    int indexTop = static_cast<int>(floor((topZ-origin[2])/spacing[2]));
	    int indexBottom = static_cast<int>(ceil((bottomZ-origin[2])/spacing[2]));

	    // set edge voxels
	    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,indexTop));
	    p[0] = boundVal;
	    ijk[2] = indexTop;
	    subList1[myThread]->InsertNextId(breast->ComputePointId(ijk));		
	    p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,indexBottom));
	    p[0] = boundVal;
	    ijk[2] = indexBottom;
	    if(indexBottom != indexTop){
	      subList1[myThread]->InsertNextId(breast->ComputePointId(ijk));
	    }
	    // set voxels between these 2 points to inner value;
	    for(int k=indexTop+1; k<indexBottom; k++){
	      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	      p[0] = innerVal;
	    }
	  }
	}
      }
    }

    // cout << "Thread " << myThread << ": direction 1 done\n";

    // find intersect with second set of directions
    // iterate over x and z values
    for(int i=myThread; i<dim[0]; i+=numThread){
      double xpos = origin[0]+i*spacing[0];
      int ijk[3];
      ijk[0] = i;
      for(int j=0; j<dim[2]; j++){

	double zpos = origin[2]+j*spacing[2];
	ijk[2] = j;
	// calculate y position of top surface
	double lineStart[3]; // end points of line
	double lineEnd[3];
      
	double tol = 0.005;
	double tval; 
	vtkIdType intersectCell;
	int subId;  

	double intersect[3]; // output position
	double pcoords[3];

	lineStart[0] = xpos;
	lineStart[1] = baseBound[2];
	lineStart[2] = zpos;

	lineEnd[0] = xpos;
	lineEnd[1] = baseBound[3];
	lineEnd[2] = zpos;

	if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					   tval, intersect, pcoords, subId, intersectCell)){

	  // found intersection
	  double topY = intersect[1];

	  // do other direction
	  lineStart[1] = baseBound[3];
	  lineEnd[1] = baseBound[2];

	  if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					     tval, intersect, pcoords, subId, intersectCell)){
	    
	    double bottomY = intersect[1];

	    // find nearest voxels to intersections
	    int indexTop = static_cast<int>(floor((topY-origin[1])/spacing[1]));
	    int indexBottom = static_cast<int>(ceil((bottomY-origin[1])/spacing[1]));
	    
	    // set edge voxels
	    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,indexTop,j));
	    p[0] = boundVal;
	    ijk[1] = indexTop;
	    subList2[myThread]->InsertNextId(breast->ComputePointId(ijk));
	    p = static_cast<unsigned char*>(breast->GetScalarPointer(i,indexBottom,j));
	    p[0] = boundVal;
	    ijk[1] = indexBottom;
	    if(indexBottom != indexTop){
	      subList2[myThread]->InsertNextId(breast->ComputePointId(ijk));
	    }
	    // set voxels between these 2 points to inner value;
	    for(int k=indexTop+1; k<indexBottom; k++){
	      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,k,j));
	      // only change if not boundary
	      if(p[0] != boundVal){
		p[0] = innerVal;
	      }
	    }
	  }
	}
      }
    }

    // cout << "Thread " << myThread << ": direction 2 done\n";

    // find intersect with final set of directions
    // iterate over y and z values
    for(int i=myThread; i<dim[1]; i+=numThread){
      double ypos = origin[1]+i*spacing[1];
      int ijk[3];
      ijk[1] = i;
      for(int j=0; j<dim[2]; j++){
      
	double zpos = origin[2]+j*spacing[2];
	ijk[2] = j;
      
	// calculate x position of top surface
	double lineStart[3]; // end points of line
	double lineEnd[3];
      
	double tol = 0.005;
	double tval; 
	vtkIdType intersectCell;
	int subId; 

	double intersect[3]; // output position
	double pcoords[3];

	lineStart[0] = baseBound[0];
	lineStart[1] = ypos;
	lineStart[2] = zpos;

	lineEnd[0] = baseBound[1];;
	lineEnd[1] = ypos;
	lineEnd[2] = zpos;

	if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					   tval, intersect, pcoords, subId, intersectCell)){
	
	  // found intersection
	  double topX = intersect[0];

	  // do other direction
	  lineStart[0] = baseBound[1];
	  lineEnd[0] = baseBound[0];

	  if(innerLocator->IntersectWithLine(lineStart, lineEnd, tol,
					     tval, intersect, pcoords, subId, intersectCell)){
	  
	    double bottomX = intersect[0];
	  
	    // find nearest voxels to intersections
	    int indexTop = static_cast<int>(floor((topX-origin[0])/spacing[0]));
	    indexTop = (indexTop < 0) ? 0 : indexTop; 
	    int indexBottom = static_cast<int>(ceil((bottomX-origin[0])/spacing[0]));
	    indexBottom = (indexBottom > dim[0]-1) ? dim[0]-1 : indexBottom;
	    indexBottom = (indexBottom < 0) ? 0 : indexBottom;

	    // set edge voxels on front side only
	    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(indexBottom,i,j));
	    p[0] = boundVal;
	    ijk[0] = indexBottom;
	    subList3[myThread]->InsertNextId(breast->ComputePointId(ijk));
	    // set voxels between these 2 points to inner value;
	    for(int k=indexTop+1; k<indexBottom; k++){
	      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(k,i,j));
	      // only change if not boundary
	      if(p[0] != boundVal){
		p[0] = innerVal;
	      }
	    }
	  }
	}
      }
    }
    // cout << "Thread " << myThread << ": direction 3 done\n";
  }
	    
  // combine lists
#pragma omp parallel sections
  {
#pragma omp section
    {
      vtkIdType nList1;
      vtkIdType c1 = 0;
      nList1 = subList1[0]->GetNumberOfIds();
      for(int i=1; i<maxThread; i++){
	nList1 += subList1[i]->GetNumberOfIds();
      }
      boundaryList1->SetNumberOfIds(nList1);
      for(int i=0; i<maxThread; i++){
	int numPts = subList1[i]->GetNumberOfIds();
	for(int j=0; j<numPts; j++){
	  boundaryList1->InsertId(c1,subList1[i]->GetId(j));
	  c1++;
	}
      }
      vtkSortDataArray::Sort(boundaryList1);
    }
#pragma omp section
    {
      vtkIdType nList2;
      vtkIdType c2 = 0;
      nList2 = subList2[0]->GetNumberOfIds();
      for(int i=1; i<maxThread;i++){
	nList2 += subList2[i]->GetNumberOfIds();
      }
      boundaryList2->SetNumberOfIds(nList2);
      for(int i=0; i<maxThread; i++){
	int numPts = subList2[i]->GetNumberOfIds();
	for(int j=0; j<numPts; j++){
	  boundaryList2->InsertId(c2,subList2[i]->GetId(j));
	  c2++;
	}
      }
      vtkSortDataArray::Sort(boundaryList2);
    }
#pragma omp section
    {
      vtkIdType nList3;
      vtkIdType c3 = 0;
      nList3 = subList3[0]->GetNumberOfIds();
      for(int i=1; i<maxThread;i++){
	nList3 += subList3[i]->GetNumberOfIds();
      }
      boundaryList3->SetNumberOfIds(nList3);
      for(int i=0; i<maxThread; i++){
	int numPts = subList3[i]->GetNumberOfIds();
	for(int j=0; j<numPts; j++){
	  boundaryList3->InsertId(c3,subList3[i]->GetId(j));
	  c3++;
	}
      }
      vtkSortDataArray::Sort(boundaryList3);
    }
  }

  vtkIdType nList1 = boundaryList1->GetNumberOfIds();
  vtkIdType nList2 = boundaryList2->GetNumberOfIds();
  vtkIdType nList3 = boundaryList3->GetNumberOfIds();

  vtkIdType* pList1 = boundaryList1->GetPointer(0);
  vtkIdType* pList2 = boundaryList2->GetPointer(0);
  vtkIdType* pList3 = boundaryList3->GetPointer(0);

  vtkIdType cList1 = 0;
  vtkIdType cList2 = 0;
  vtkIdType cList3 = 0;

  bool dList1 = false;
  bool dList2 = false;
  bool dList3 = false;

  while(!dList1 || !dList2 || !dList3){
    vtkIdType curMin;
    if(!dList1){
      curMin = *pList1;
      if(!dList2){
	curMin = *pList2 < curMin ? *pList2 : curMin;
	if(!dList3){
	  curMin = *pList3 < curMin ? *pList3 : curMin;
	  // check 1,2,3
	  boundaryList->InsertNextId(curMin);
	  while(*pList1 == curMin && !dList1){
	    cList1++;
	    if(cList1 < nList1){
	      pList1++;
	    } else {
	      dList1 = true;
	    }
	  }
	  while(*pList2 == curMin && !dList2){
	    cList2++;
	    if(cList2 < nList2){
	      pList2++;
	    } else {
	      dList2 = true;
	    }
	  }
	  while(*pList3 == curMin && !dList3){
            cList3++;
            if(cList3 < nList3){
	      pList3++;
            } else {
	      dList3 = true;
            }
          }
	} else {
	  // check 1,2
	  boundaryList->InsertNextId(curMin);
	  while(*pList1 == curMin && !dList1){
            cList1++;
            if(cList1 <nList1){
              pList1++;
	    } else {
	      dList1 = true;
            }
          }
	  while(*pList2 == curMin && !dList2){
            cList2++;
            if(cList2 <nList2){
	      pList2++;
            } else {
	      dList2 = true;
            }
          }
	}
      } else {
	if(!dList3){
	  curMin = *pList3 < curMin ? *pList3 : curMin;
	  // check 1,3
	  boundaryList->InsertNextId(curMin);
	  while(*pList1 == curMin && !dList1){
	    cList1++;
	    if(cList1 <nList1){
	      pList1++;
	    } else {
	      dList1 = true;
	    }
	  }
	  while(*pList3 == curMin && !dList3){
            cList3++;
            if(cList3 < nList3){
              pList3++;
            } else {
              dList3 = true;
            }
          }
	} else {
	  // check 1
	  boundaryList->InsertNextId(curMin);
	  while(*pList1 == curMin && !dList1){
            cList1++;
            if(cList1 <nList1){
              pList1++;
            } else {
              dList1 = true;
            }
          }
	}
      }
    } else {
      // 1 done
      if(!dList2){
	curMin = *pList2;
	if(!dList3){
	  curMin = *pList3 < curMin ? *pList3 : curMin;
	  // check 2,3
	  boundaryList->InsertNextId(curMin);
	  while(*pList2 == curMin && !dList2){
            cList2++;
            if(cList2 <nList2){
              pList2++;
            } else {
              dList2 = true;
            }
          }
	  while(*pList3 == curMin && !dList3){
            cList3++;
            if(cList3 < nList3){
              pList3++;
            } else {
              dList3 = true;
            }
          }
	} else {
	  // check 2
	  boundaryList->InsertNextId(curMin);
	  while(*pList2 == curMin && !dList2){
            cList2++;
            if(cList2 <nList2){
              pList2++;
            } else {
              dList2 = true;
            }
          }
	}
      } else {
	// only 3 left
	curMin = *pList3;
	boundaryList->InsertNextId(curMin);
	while(*pList3 == curMin && !dList3){
	  cList3++;
	  if(cList3 < nList3){
	    pList3++;
	  } else {
	    dList3 = true;
	  }
	}
      }
    }
  }
	    
  // delete array subList
  delete [] subList1;
  delete [] subList2;
  delete [] subList3;

  //cout << "done.\n";

  // correct boundary list to be all boundary values 
  vtkIdType dnum = boundaryList->GetNumberOfIds();
  for(vtkIdType i=0; i<dnum; i++){
    double loc[3];
    int ijk[3];
    double pcoords[3];
    breast->GetPoint(boundaryList->GetId(i),loc);
    breast->ComputeStructuredCoordinates(loc,ijk,pcoords);
    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(ijk));
    if(p[0] != boundVal){
      p[0] = boundVal;
    }
  }

  /***********************
	Skin
  ***********************/
  // timing
  //time(&currentTime);
  //cout << "Time interval 2: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  // grow a skin layer
  //cout << "Growing " << skinThick << " mm skin layer...";

  // create list of surrounding voxels to check
  vtkSmartPointer<vtkIntArray> checkVoxels =
    vtkSmartPointer<vtkIntArray>::New();

  checkVoxels->SetNumberOfComponents(3);

  int voxelThick = (int)ceil(skinThick/imgRes);

  for(int i=-1*voxelThick; i<=voxelThick; i++){
    for(int j=-1*voxelThick; j<=voxelThick; j++){
      for(int k=-1*voxelThick; k<=voxelThick; k++){
	if(sqrt((double)(i*i+j*j+k*k))*imgRes <= skinThick){
	  // found a voxel to check
	  checkVoxels->InsertNextTuple3(i,j,k);
	}
      }
    }
  }

  vtkIdType numCheck = checkVoxels->GetNumberOfTuples();

  int breastExtent[6];
  breast->GetExtent(breastExtent);

  // breast volume in voxels
  long long int breastVoxVol = 0;

  // only add skin for x>0
  int minSkinXVox = static_cast<int>(ceil(-origin[0]/imgRes));
  int areolaVoxel = static_cast<int>(ceil(areolaRad/imgRes));

  // add skin thickness near nipple
  double skinThick2 = skinThick*2.0;

  // find point on mesh closest to center
  vtkSmartPointer<vtkPointLocator> locator =
    vtkSmartPointer<vtkPointLocator>::New();
  locator->SetDataSet(innerPoly);
  locator->BuildLocator();

  vtkIdType nipplePt;
  nipplePt = locator->FindClosestPoint(nipplePos);

  int nippleVoxel[3];   // coordinates of nipple base
  double nipplePCoords[3]; // parametric coordinates
  breast->ComputeStructuredCoordinates(nipplePos,nippleVoxel,nipplePCoords);

  // iterate over boundary voxel list, grow skin
  vtkIdType nCurBoundary = boundaryList->GetNumberOfIds();

  //vtkSmartPointer<vtkIdList> boundaryDel2 =
  //  vtkSmartPointer<vtkIdList>::New();

  maxThread = omp_get_max_threads();

  // boundary sub-list
  //vtkSmartPointer<vtkIdList> *boundaryDelSub2 = new vtkSmartPointer<vtkIdList>[maxThread];
  //for(int i=0; i<maxThread; i++){
  //  boundaryDelSub2[i] = vtkSmartPointer<vtkIdList>::New();
  //}

#pragma omp parallel for num_threads(maxThread)
  for(vtkIdType i=0; i<nCurBoundary; i++){
    int myThread = omp_get_thread_num();
    vtkIdType myId = boundaryList->GetId(i);
    double loc[3];
    double pcoords[3];
    int ijk[3];
    breast->GetPoint(myId, loc);
    breast->ComputeStructuredCoordinates(loc, ijk, pcoords);

    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(ijk));
    
    //if(p[0] == boundVal){
      // good boundary point
    if(ijk[0] >= minSkinXVox){
      double nipDist2 = vtkMath::Distance2BetweenPoints(loc,nipplePos);
      if(nipDist2 > 4*areolaRad*areolaRad){
	// boundary voxel for skinning
	for(vtkIdType m=0; m<numCheck; m++){
	  double offset[3];
	  checkVoxels->GetTuple(m,offset);
	  int a,b,c;
	  a = ijk[0]+(int)offset[0];
	  b = ijk[1]+(int)offset[1];
	  c = ijk[2]+(int)offset[2];
	  if(a>=breastExtent[0] && a<=breastExtent[1] && b>=breastExtent[2] && b<=breastExtent[3] &&
	     c>=breastExtent[4] && c<=breastExtent[5]){
	    unsigned char* q =
	      static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
	    if(q[0] == tissue.bg){
#pragma omp atomic write
	      q[0] = tissue.skin;
	    }
	  }
	}
      } else {
	// areola
	double mySkinThick = skinThick + (skinThick2-skinThick)/(1+exp(12/areolaRad*(sqrt(nipDist2)-areolaRad)));
	int mySearchRad = static_cast<int>(ceil(mySkinThick/imgRes));
	for(int a=ijk[0]-mySearchRad; a<=ijk[0]+mySearchRad; a++){
	  for(int b=ijk[1]-mySearchRad; b<=ijk[1]+mySearchRad; b++){
	    for(int c=ijk[2]-mySearchRad; c<=ijk[2]+mySearchRad; c++){
	      unsigned char* q = static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
	      if(q[0] == tissue.bg){
		// check distance                                                                                                                                                                                                
		double skinDist = imgRes*sqrt(static_cast<double>((a-ijk[0])*(a-ijk[0])+(b-ijk[1])*(b-ijk[1])+(c-ijk[2])*(c-ijk[2])));
		if(skinDist <= mySkinThick){
#pragma omp atomic write
		  q[0] = tissue.skin;
		}
	      }
	    }
	  }
	}
      }
    }
    p[0] = innerVal;
    //} else {
    //  // point should not be in boundary
    //  boundaryDelSub2[myThread]->InsertNextId(myId);
    //}
  }

  // update boundary
  //for(int i=0; i<maxThread; i++){
  //  vtkIdType nDel = boundaryDelSub2[i]->GetNumberOfIds();
  //  // debug
  //  cout << "Delete list" << i << " length = " << nDel << "\n";
  //  for(int j=0; j<nDel; j++){
  //    boundaryList->DeleteId(boundaryDelSub2[i]->GetId(j));
  //  }
  //}

  // calculate inner volume and correct border errors
#pragma omp parallel for reduction(+:breastVoxVol)
  for(int i=0; i<dim[0]; i++){
    int ijk[3];
    ijk[0] = i;
    for(int j=0; j<dim[1]; j++){
      ijk[1] = j;
      for(int k=0; k<dim[2]; k++){
	ijk[2] = k;
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(ijk));
	if(p[0] == innerVal){
	  breastVoxVol += 1;
	} else {
	  if(p[0] == boundVal){
	    p[0] = innerVal;
	    breastVoxVol += 1;
	  }
	}
      }
    }
  }

  /*
  // add skin thickness near nipple
  double skinThick2 = skinThick*2.0;
	
  // find point on mesh closest to center
  vtkSmartPointer<vtkPointLocator> locator =
    vtkSmartPointer<vtkPointLocator>::New();
  locator->SetDataSet(innerPoly);
  locator->BuildLocator();

  vtkIdType nipplePt;
  nipplePt = locator->FindClosestPoint(nipplePos);
	
  int nippleVoxel[3];	// coordinates of nipple base
  double nipplePCoords[3]; // parametric coordinates
  breast->ComputeStructuredCoordinates(nipplePos,nippleVoxel,nipplePCoords);
	
  for(int i=nippleVoxel[0]-2*areolaVoxel; i<=dim[0]-1; i++){
    for(int j=nippleVoxel[1]-2*areolaVoxel; j<=nippleVoxel[1]+2*areolaVoxel; j++){
      for(int k=nippleVoxel[2]-2*areolaVoxel; k<=nippleVoxel[2]+2*areolaVoxel; k++){
	int ijk[3] = {i,j,k};
	double loc[3];
	breast->GetPoint(breast->ComputePointId(ijk),loc);
	double nipDist2 = vtkMath::Distance2BetweenPoints(loc,nipplePos);
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == boundVal && nipDist2 <= 4*areolaRad*areolaRad){
	  double mySkinThick = skinThick + (skinThick2-skinThick)/(1+exp(12/areolaRad*(sqrt(nipDist2)-areolaRad)));
	  int mySearchRad = static_cast<int>(ceil(mySkinThick/imgRes));
	  for(int a=i-mySearchRad; a<=i+mySearchRad; a++){
	    for(int b=j-mySearchRad; b<=j+mySearchRad; b++){
	      for(int c=k-mySearchRad; c<=k+mySearchRad; c++){
		unsigned char* q = static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
		if(q[0] == tissue.bg){
		  // check distance
		  double skinDist = imgRes*sqrt(static_cast<double>((a-i)*(a-i)+(b-j)*(b-j)+(c-k)*(c-k)));
		  if(skinDist <= mySkinThick){
		    q[0] = tissue.skin;
		  }
		}
	      }
	    }
	  }
	  // can now change boundary to interior
	  p[0] = innerVal;
	  breastVoxVol += 1;
	}
      }
    }
  }					

  */
				
  double breastVol = (double)breastVoxVol*pow(imgRes,3.0);

  //cout << "done.\n";

  //cout << "Breast volume: " << breastVol/1000 << " cc ("<< breastVoxVol << " voxels).\n";

  /***********************
	Nipple
  **********************/

  // debug
  //time(&currentTime);
  //cout << "Time interval 9: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  // create nipple structure
  //cout << "Creating nipple structure...";

  // squared radius
  double nippleRad2 = nippleRad*nippleRad;
	
  double nippleNorm[3];
  // extract normal
  vtkSmartPointer<vtkFloatArray> normals2 =
    vtkFloatArray::SafeDownCast(innerPoly->GetPointData()->GetNormals());

  normals2->GetTuple(nipplePt, nippleNorm);
  double nippleNormLen = sqrt(nippleNorm[0]*nippleNorm[0] +
			      nippleNorm[1]*nippleNorm[1] + nippleNorm[2]*nippleNorm[2]);
  for(int i=0; i<3; i++){
    nippleNorm[i] = nippleNorm[i]/nippleNormLen;
  }
  // nippleNorm may be inward pointing
  if(nippleNorm[0] < 0.0){
    for(int i=0; i<3; i++){
      nippleNorm[i] = -1.0*nippleNorm[i];
    }
  }

  // create nipple

  // nipple function
  // superquadric (rad/nippleRad)^t+abs(len/nippleLen)^t <= 1 t = 2.5 - 8
  double nippleShape = 3.0;
  
  double searchRad = sqrt(2*nippleRad*nippleRad+nippleLen*nippleLen);

  int searchSpace[6];	// search region for nipple
  searchSpace[0] = nippleVoxel[0] - (int)ceil(searchRad/spacing[0]);
  if(searchSpace[0] < breastExtent[0]){
    searchSpace[0] = breastExtent[0];
  }
  searchSpace[1] = nippleVoxel[0] + (int)ceil(searchRad/spacing[0]);
  if(searchSpace[1] > breastExtent[1]){
    searchSpace[1] = breastExtent[1];
  }
  searchSpace[2] = nippleVoxel[1] - (int)ceil(searchRad/spacing[1]);
  if(searchSpace[2] < breastExtent[2]){
    searchSpace[2] = breastExtent[2];
  }
  searchSpace[3] = nippleVoxel[1] + (int)ceil(searchRad/spacing[1]);
  if(searchSpace[3] > breastExtent[3]){
    searchSpace[3] = breastExtent[3];
  }
  searchSpace[4] = nippleVoxel[2] - (int)ceil(searchRad/spacing[2]);
  if(searchSpace[4] < breastExtent[4]){
    searchSpace[4] = breastExtent[4];
  }
  searchSpace[5] = nippleVoxel[2] + (int)ceil(searchRad/spacing[2]);
  if(searchSpace[5] > breastExtent[5]){
    searchSpace[5] = breastExtent[5];
  }

#pragma omp parallel for
  for(int i=searchSpace[0]; i<= searchSpace[1]; i++){
    for(int j=searchSpace[2]; j<= searchSpace[3]; j++){
      for(int k=searchSpace[4]; k<= searchSpace[5]; k++){
	// search cube, project onto normal and check if in nipple bound
	// only if outside breast interior
	unsigned char* q =
	  static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(q[0] != innerVal){ // not in interior
	  // get id for point
	  vtkIdType id;
	  int coord[3];
	  coord[0] = i;
	  coord[1] = j;
	  coord[2] = k;
	  id = breast->ComputePointId(coord);
	  // get spatial coordinates of point
	  double pos[3];
	  breast->GetPoint(id,pos);
	  // compute distance to nipple line and length along line
	  double dist=0.0;
	  double len=0.0;
	  // project onto nipple line
	  for(int m=0; m<3; m++){
	    len += nippleNorm[m]*(pos[m]-nipplePos[m]);
	  }
	  // distance from nipple line
	  for(int m=0; m<3; m++){
	    dist += (pos[m]-nipplePos[m]-len*nippleNorm[m])*(pos[m]-nipplePos[m]-len*nippleNorm[m]);
	  }
	  dist = sqrt(dist);
	  
	  if(pow(dist/nippleRad,nippleShape)+pow(fabs(len)/nippleLen,nippleShape) <= 1.0){
	    q[0] = tissue.nipple;
	  }
	}
      }
    }
  }

  //cout << "done.\n";

  // add chest muscle

#pragma omp parallel for  
  for(int j=0; j<dim[1]; j++){
	
    int muscleThick;
		
    if(leftSide){
      muscleThick = static_cast<int>(ceil((minSkinXVox-1)*(1-static_cast<double>(j*j)/(dim[1]*dim[1]))));
    }else{
      muscleThick = static_cast<int>(ceil((minSkinXVox-1)*(1-static_cast<double>((j-dim[1])*(j-dim[1]))/(dim[1]*dim[1]))));
    }
		
    for(int k=0; k<dim[2]; k++){
      for(int i=0; i<=muscleThick; i++){
	unsigned char* p =
	  static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == innerVal){
	  p[0] = tissue.muscle;
	}
      }
    }
  }

  // save the volumetric data
  // save volume image data
  //vtkSmartPointer<vtkXMLImageDataWriter> writer2 =
  //  vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //writer2->SetFileName(breastVoxelFilename.c_str());
	
  //#if VTK_MAJOR_VERSION <= 5
  //writer2->SetInput(breast);
  //#else
  //writer2->SetInputData(breast);
  //#endif
  //writer2->Write();

  /***********************
	Compartments
  ***********************/

  // timing
  //time(&currentTime);
  //cout << "Time interval 3: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  // breast segmentation into compartments and lipid buffer zone
  //cout << "Generating breast compartment seed points...";

  int numBreastCompartments = vm["compartments.num"].as<int>();
  
  int numAngles = numBreastCompartments;
	
  double seedBaseDist = vm["compartments.seedBaseDist"].as<double>();	// distance along nipple line of seed base

  //double backFatBufferFrac = vm["compartments.backFatBufferFrac"].as<double>(); // fraction of phantom in x-direction forced to be fat
  // pick backplane seed points on this plane
  int numBackSeeds = vm["compartments.numBackSeeds"].as<int>();	// number of backplane seed points
  int numSkinSeeds = 250;

  double angularJitter = 2*pi*numAngles*vm["compartments.angularJitter"].as<double>();	// max angular jitter for seed placement
  double zJitter = vm["compartments.zJitter"].as<double>();	// jitter in z direction
  double maxFracRadialDist = vm["compartments.maxFracRadialDist"].as<double>(); // minimum and maximum radial distance from baseSeed as fraction of distance to breast surface
  double minFracRadialDist = vm["compartments.minFracRadialDist"].as<double>();

  double scaleMin[3];  // scaling for gland compartments
  scaleMin[0] = vm["compartments.minScaleNippleDir"].as<double>();
  scaleMin[1] = vm["compartments.minScale"].as<double>();
  scaleMin[2] = vm["compartments.minScale"].as<double>();
	
  double scaleMax[3];	// 1st dimension principally points toward nipple - gland seeds
  scaleMax[0] = vm["compartments.maxScaleNippleDir"].as<double>();
  scaleMax[1] = vm["compartments.maxScale"].as<double>();
  scaleMax[2] = vm["compartments.maxScale"].as<double>();
	
  double gMin = vm["compartments.minGlandStrength"].as<double>();	// strength of gland compartments
  double gMax = vm["compartments.maxGlandStrength"].as<double>();
  double deflectMax = pi*vm["compartments.maxDeflect"].as<double>();	// maximum deflection angle from pointing towards nipple (gland compartments)

  // backplane, nipple are spherical weighting
  // skin orients towards nipple
  double scaleSkinMin[3];
  scaleSkinMin[0] = vm["compartments.minSkinScaleNippleDir"].as<double>();
  scaleSkinMin[1] = vm["compartments.minSkinScale"].as<double>();
  scaleSkinMin[2] = vm["compartments.minSkinScale"].as<double>();
  double scaleSkinMax[3];
  scaleSkinMax[0] = vm["compartments.maxSkinScaleNippleDir"].as<double>();
  scaleSkinMax[1] = vm["compartments.maxSkinScale"].as<double>();
  scaleSkinMax[2] = vm["compartments.maxSkinScale"].as<double>();
	
  double gSkin = vm["compartments.skinStrength"].as<double>();

  double scaleBack = vm["compartments.backScale"].as<double>();
  double gBack = vm["compartments.backStrength"].as<double>();

  double scaleNipple = vm["compartments.nippleScale"].as<double>();
  double gNipple = vm["compartments.nippleStrength"].as<double>();
	
  // debug
  double ligDistThresh = 12.0;
	
  // radius of fat seeds to check for segmentation
  double compSeedRadius = vm["compartments.voronSeedRadius"].as<double>();

  int numFatSeeds = numAngles + numBackSeeds + numSkinSeeds;  // total number of fat seed points

  typedef struct breastComp{
    // data structure for breast compartment
    double scale[3]; // size
    double g;	// size weight
    vtkVector3d axis[3];	// local coordinate system
    double pos[3];	// position
    bool fat;	// is it fat?
    bool keep;  // retain glandular compartment?
    unsigned char compId;	// compartment id 1,2,....
    int boundBox[6];	// bounding box indices
    int voxelCount;		// number of voxels
    double volume;		// volume (cubic mm)
  } breastComp;

  // memory allocation
  breastComp *glandCompartments = (breastComp*)malloc((numBreastCompartments+1)*sizeof(breastComp));
  breastComp *fatCompartments = (breastComp*)malloc(numFatSeeds*sizeof(breastComp));
	
  // breast dimensions for initializing bounding box
  int breastDim[3];
  breast->GetDimensions(breastDim);

  vtkSmartPointer<vtkPoints> seeds =
    vtkSmartPointer<vtkPoints>::New();

  // base of spokes for seed placement
  double seedBase[3];
  // coordinate system for seedBase
  vtkVector3d baseAxis[3];

  // set seed base location and specify first coordinate vector
  for(int i=0; i<3; i++){
    seedBase[i] = nipplePos[i]-seedBaseDist*nippleNorm[i];
    baseAxis[0][i] = -1.0*nippleNorm[i];
  }

  // check seed base within breast
  int coords[3];	// coordinates of nipple seed base
  double pcoords[3]; // parametric coordinates
  breast->ComputeStructuredCoordinates(seedBase,coords,pcoords);
  unsigned char* base = static_cast<unsigned char*>(breast->GetScalarPointer(coords));
  if(base[0] != innerVal){
    // outside breast error
    cout << "Error, breast compartment seed base outside breast volume\n";
    return EXIT_FAILURE;
  }

  // construct other coordinate vectors

  // second coordinate vector based on Gram-Schmidt using (0,1,0)
  vtkVector3d v2;
  v2[0] = seedBase[0];
  v2[1] = seedBase[1] - 1.0;
  v2[2] = seedBase[2];
  double innerProd = v2.Dot(baseAxis[0]);

  for(int i=0; i<3; i++){
    baseAxis[1][i] = v2[i] - innerProd*baseAxis[0][i];
  }
  // normalize
  baseAxis[1].Normalize();

  // calculate 3rd vector based on cross product
  baseAxis[2] = baseAxis[0].Cross(baseAxis[1]);

  /* now have baseAxis coordinate system
     shoot rays toward breast surface and determine compartment seed locations
     and breast surface seed points */
  
  double rayLength = 500.0;	// some large value to guarantee being outside of breast

  for(int i=0; i<numAngles; i++){
    double theta = (double)i*2*pi/(double)numAngles;
    
    // add theta jitter
    theta = theta + rgen->GetRangeValue(-1.0*angularJitter,angularJitter);
    rgen->Next();

    vtkVector3d rayDir;
    double lineEnd[3];
    for(int j=0; j<3; j++){
      rayDir[j] = cos(theta)*baseAxis[1][j]+sin(theta)*baseAxis[2][j];
      lineEnd[j] = seedBase[j]+rayLength*rayDir[j];
    }

    // intersection variables
    double tol = 0.005;
    double tval; // not sure what this is for
    vtkIdType intersectCell;
    int subId;  // probably don't need this

    double intersect[3]; // output position
    double pcoords[3];

    vtkSmartPointer<vtkCellLocator> innerLocator =
      vtkSmartPointer<vtkCellLocator>::New();
    innerLocator->SetDataSet(innerPoly);
    innerLocator->BuildLocator();

    if(innerLocator->IntersectWithLine(seedBase, lineEnd, tol,
				       tval, intersect, pcoords, subId, intersectCell)){
      
      // the breast surface seed point is intersect
      seeds->InsertNextPoint(intersect);
      
      // add fat seed to structure
      for(int j=0; j<3; j++){
	fatCompartments[i].pos[j] = intersect[j];
	fatCompartments[i].scale[j] = rgen->GetRangeValue(scaleSkinMin[j],scaleSkinMax[j]);
	rgen->Next();
      }
      fatCompartments[i].g = gSkin;
      fatCompartments[i].fat = true;
      fatCompartments[i].keep = true;
      fatCompartments[i].compId = 0;
      fatCompartments[i].voxelCount = 0;

      fatCompartments[i].boundBox[0] = breastDim[0]+1;
      fatCompartments[i].boundBox[1] = -1;
      fatCompartments[i].boundBox[2] = breastDim[1]+1;
      fatCompartments[i].boundBox[3] = -1;
      fatCompartments[i].boundBox[4] = breastDim[2]+1;
      fatCompartments[i].boundBox[5] = -1;

      // principle unit vector oriented towards nipple
      // first vector
      for(int j=0; j<3; j++){
	fatCompartments[i].axis[0][j] = nipplePos[j] - intersect[j];
      }
      // normalize
      fatCompartments[i].axis[0].Normalize();

      // calculate second vector based on direction to coordinate origin
      vtkVector3d v2;
      for(int j=0; j<3; j++){
	v2[j] = intersect[j];
      }
      double innerProd = v2.Dot(fatCompartments[i].axis[0]);

      for(int j=0; j<3; j++){
	fatCompartments[i].axis[1][j] = v2[j] - innerProd*fatCompartments[i].axis[0][j];
      }
      // normalize
      fatCompartments[i].axis[1].Normalize();
			
      // calculate 3rd vector based on cross product
      fatCompartments[i].axis[2] = fatCompartments[i].axis[0].Cross(fatCompartments[i].axis[1]);

      // find gland compartment seed point

      double skinDist = sqrt(vtkMath::Distance2BetweenPoints(seedBase,intersect));
      double seedDist = skinDist*rgen->GetRangeValue(minFracRadialDist,
						     maxFracRadialDist);
      rgen->Next();

      double seedPos[3];
      for(int j=0; j<3; j++){
	seedPos[j] = seedBase[j]+seedDist*rayDir[j];
      }

      // add z jitter
      
      double zjit = rgen->GetRangeValue(-1.0*zJitter, 1.0*zJitter);
      rgen->Next();

      for(int j=0; j<3; j++){
	seedPos[j] = seedPos[j]+zjit*baseAxis[0][j];
      }

      //seeds->InsertNextPoint(seedPos);

      // add to gland structure
      for(int j=0; j<3; j++){
	glandCompartments[i].pos[j] = seedPos[j];
	glandCompartments[i].scale[j] = rgen->GetRangeValue(scaleMin[j],scaleMax[j]);
	rgen->Next();
      }
			
      glandCompartments[i].g = rgen->GetRangeValue(gMin,gMax);
      rgen->Next();
      glandCompartments[i].fat = false;
      glandCompartments[i].keep = true;
      glandCompartments[i].compId = (unsigned char)i;
      glandCompartments[i].voxelCount = 0;

      glandCompartments[i].boundBox[0] = breastDim[0]+1;
      glandCompartments[i].boundBox[1] = -1;
      glandCompartments[i].boundBox[2] = breastDim[1]+1;
      glandCompartments[i].boundBox[3] = -1;
      glandCompartments[i].boundBox[4] = breastDim[2]+1;
      glandCompartments[i].boundBox[5] = -1;

      // coordinate system
      // first vector
      for(int j=0; j<3; j++){
	glandCompartments[i].axis[0][j] = nipplePos[j] - seedPos[j];
      }
      // normalize
      glandCompartments[i].axis[0].Normalize();

      // calculate second vector based on direction to coordinate origin
      for(int j=0; j<3; j++){
	v2[j] = seedPos[j];
      }
      innerProd = v2.Dot(glandCompartments[i].axis[0]);
      
      for(int j=0; j<3; j++){
	glandCompartments[i].axis[1][j] = v2[j] - innerProd*glandCompartments[i].axis[0][j];
      }
      // normalize
      glandCompartments[i].axis[1].Normalize();

      // calculate 3rd vector based on cross product
      glandCompartments[i].axis[2] = glandCompartments[i].axis[0].Cross(glandCompartments[i].axis[1]);

      // have 3 unit vectors
      // rotate randomly about principle direction (to nipple)
      double dtheta = rgen->GetRangeValue(0,2*pi);
      rgen->Next();
      double dphi = rgen->GetRangeValue(0,deflectMax);
      rgen->Next();
      double dr = tan(pi*dphi);

      for(int j=0; j<3; j++){
	glandCompartments[i].axis[0][j] = nipplePos[j] - seedPos[j] + dr*cos(dtheta)*glandCompartments[i].axis[1][j] +
	  dr*sin(dtheta)*glandCompartments[i].axis[2][j];
      }
      // normalize
      glandCompartments[i].axis[0].Normalize();
      
      // re-calculate second vector based on direction to coordinate origin and updated principle direction
      for(int j=0; j<3; j++){
	v2[j] = seedPos[j];
      }
      innerProd = v2.Dot(glandCompartments[i].axis[0]);
      
      for(int j=0; j<3; j++){
	glandCompartments[i].axis[1][j] = v2[j] - innerProd*glandCompartments[i].axis[0][j];
      }
      // normalize
      glandCompartments[i].axis[1].Normalize();

      // calculate 3rd vector based on cross product
      glandCompartments[i].axis[2] = glandCompartments[i].axis[0].Cross(glandCompartments[i].axis[1]);

      // finished calculating axes
    } else {
      // missed the breast surface - error
      cout << "Error, missed breast surface when shooting rays\n";
      return EXIT_FAILURE;
    }
  }

  // add nipple
  // seeds->InsertNextPoint(nipplePos);
  
  // debug
  double nipSc[3] = {20, 1.0, 1.0};

  // add nipple seed to gland structure
  for(int j=0; j<3; j++){
    glandCompartments[numBreastCompartments].pos[j] = nipplePos[j];
    glandCompartments[numBreastCompartments].scale[j] = scaleNipple;
  }
  glandCompartments[numBreastCompartments].g = gNipple;
  glandCompartments[numBreastCompartments].fat = false;
  glandCompartments[numBreastCompartments].keep = true;
  glandCompartments[numBreastCompartments].compId = (unsigned char)numBreastCompartments;
  glandCompartments[numBreastCompartments].voxelCount = 0;

  glandCompartments[numBreastCompartments].boundBox[0] = breastDim[0]+1;
  glandCompartments[numBreastCompartments].boundBox[1] = -1;
  glandCompartments[numBreastCompartments].boundBox[2] = breastDim[1]+1;
  glandCompartments[numBreastCompartments].boundBox[3] = -1;
  glandCompartments[numBreastCompartments].boundBox[4] = breastDim[2]+1;
  glandCompartments[numBreastCompartments].boundBox[5] = -1;

  // first axis is nipple direction
  for(int k=0; k<3; k++){
    glandCompartments[numBreastCompartments].axis[0][k] = nippleNorm[k];
  }

  // others are random
  vtkMath::Perpendiculars(glandCompartments[numBreastCompartments].axis[0].GetData(),
			  glandCompartments[numBreastCompartments].axis[1].GetData(), glandCompartments[numBreastCompartments].axis[2].GetData(),
			  rgen->GetRangeValue(0, 2*pi));
  rgen->Next();
	
  // add backplane seed points
  // random points on backplane
  // backplane just above muscle layer
  double backPlanePos = 0.0 + spacing[0];
  int backPlaneInd = static_cast<int>(ceil((backPlanePos-origin[0])/spacing[0]));

  for(int i=0; i<numBackSeeds; i++){
    bool foundSeed = false;
    double y,z;
    while(!foundSeed){
      // pick random location on backplane away from edge of voxel space
      y = rgen->GetRangeValue(baseBound[2]+2*spacing[1],baseBound[3]-2*spacing[1]);
      rgen->Next();
      z = rgen->GetRangeValue(baseBound[4]+2*spacing[2],baseBound[5]-2*spacing[2]);
      rgen->Next();

      // find nearest voxel and test if in breast interior
      int yInd = static_cast<int>(floor((y-origin[1])/spacing[1]));
      int zInd = static_cast<int>(floor((z-origin[2])/spacing[2]));

      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(backPlaneInd,yInd,zInd));
      if(p[0] == innerVal){
	// found a new seed point
	foundSeed = true;
      }
    }

    // add new seed point
    double backSeed[3];
    backSeed[0] = backPlanePos;
    backSeed[1] = y;
    backSeed[2] = z;

    seeds->InsertNextPoint(backSeed);

    // add to structure
    for(int j=0; j<3; j++){
      fatCompartments[numAngles+i].pos[j] = backSeed[j];
      fatCompartments[numAngles+i].scale[j] = scaleBack;
    }
    fatCompartments[numAngles+i].g = gBack;
    fatCompartments[numAngles+i].fat = true;
    fatCompartments[numAngles+i].keep = true;
    fatCompartments[numAngles+i].compId = 0;
    fatCompartments[numAngles+i].voxelCount = 0;

    fatCompartments[numAngles+i].boundBox[0] = breastDim[0]+1;
    fatCompartments[numAngles+i].boundBox[1] = -1;
    fatCompartments[numAngles+i].boundBox[2] = breastDim[1]+1;
    fatCompartments[numAngles+i].boundBox[3] = -1;
    fatCompartments[numAngles+i].boundBox[4] = breastDim[2]+1;
    fatCompartments[numAngles+i].boundBox[5] = -1;

    // unit vectors are standard
    for(int j=0; j<3; j++){
      for(int k=0; k<3; k++){
	fatCompartments[numAngles+i].axis[j][k] = 0.0;
      }
    }
    fatCompartments[numAngles+i].axis[0][0] = 1.0;
    fatCompartments[numAngles+i].axis[1][1] = 1.0;
    fatCompartments[numAngles+i].axis[2][2] = 1.0;
  }
  // finished adding backplane seed points

  // random skin seed points
  for(int i=0; i<numSkinSeeds; i++){
    bool foundSeed = false;
    vtkIdType numPts = innerPoly->GetNumberOfPoints();
    double seedCoords[3];
    while(!foundSeed){
      // pick random breast surface point
      vtkIdType tryId = static_cast<vtkIdType>(ceil(rgen->GetRangeValue(0, numPts-1)));
      rgen->Next();
      innerPoly->GetPoint(tryId, seedCoords);
			
      // select point if not to close to backplane or nipple
      if(seedCoords[0] > backPlanePos + 10.0 && seedCoords[0] < nipplePos[0] - 10.0){
	foundSeed = true;
      }  
    }

    seeds->InsertNextPoint(seedCoords);
    // debug
    //double skinSc[3] = {150.0, 100.0, 100.0}; 

    // add to structure
    for(int j=0; j<3; j++){
      fatCompartments[numAngles+numBackSeeds+i].pos[j] = seedCoords[j];
      // debug
      //fatCompartments[numAngles+numBackSeeds+i].scale[j] = skinSc[j];
      fatCompartments[numAngles+numBackSeeds+i].scale[j] = rgen->GetRangeValue(scaleSkinMin[j],scaleSkinMax[j]);
      rgen->Next();
    }
    fatCompartments[numAngles+numBackSeeds+i].g = gSkin;
    fatCompartments[numAngles+numBackSeeds+i].fat = true;
    fatCompartments[numAngles+numBackSeeds+i].keep = true;
    fatCompartments[numAngles+numBackSeeds+i].compId = 0;
    fatCompartments[numAngles+numBackSeeds+i].voxelCount = 0;
    
    fatCompartments[numAngles+numBackSeeds+i].boundBox[0] = breastDim[0]+1;
    fatCompartments[numAngles+numBackSeeds+i].boundBox[1] = -1;
    fatCompartments[numAngles+numBackSeeds+i].boundBox[2] = breastDim[1]+1;
    fatCompartments[numAngles+numBackSeeds+i].boundBox[3] = -1;
    fatCompartments[numAngles+numBackSeeds+i].boundBox[4] = breastDim[2]+1;
    fatCompartments[numAngles+numBackSeeds+i].boundBox[5] = -1;

    // principle unit vector normal to skin
    normals2->GetTuple(nipplePt, fatCompartments[numAngles+numBackSeeds+i].axis[0].GetData());
    fatCompartments[numAngles+numBackSeeds+i].axis[0].Normalize();

    // other directions random
    vtkMath::Perpendiculars(fatCompartments[numAngles+numBackSeeds+i].axis[0].GetData(),
			    fatCompartments[numAngles+numBackSeeds+i].axis[1].GetData(), fatCompartments[numAngles+numBackSeeds+i].axis[2].GetData(),
			    rgen->GetRangeValue(0, 2*pi));
    rgen->Next();
  }

  // debug
  // create Perlin noise distance function for glandular compartments
  perlinNoise *boundary = static_cast<perlinNoise*>(::operator new(sizeof(perlinNoise)*(numBreastCompartments+1)));

  for(int i=0; i<=numBreastCompartments; i++){
    new(&boundary[i]) perlinNoise(vm, (int32_t)rgen->GetRangeValue(-1073741824, 1073741824),"boundary");
    rgen->Next();
  }

  //cout << "done.\n";

  //cout << "Segmenting glandular compartments...";

  // setup point locator
  vtkSmartPointer<vtkOctreePointLocator> findSeed =
    vtkSmartPointer<vtkOctreePointLocator>::New();
	
  // create a vtkPolyData
  vtkSmartPointer<vtkPolyData> seedSet =
    vtkSmartPointer<vtkPolyData>::New();
  seedSet->SetPoints(seeds);

  // add verticies
  vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter2 =
    vtkSmartPointer<vtkVertexGlyphFilter>::New();
#if VTK_MAJOR_VERSION <= 5
  vertexGlyphFilter2->AddInput(seedSet);
#else
  vertexGlyphFilter2->AddInputData(seedSet);
#endif
  vertexGlyphFilter2->Update();

  findSeed->SetDataSet(vertexGlyphFilter2->GetOutput());
  findSeed->BuildLocator();

  // iterate over voxels to do segmentation
	
  // starting by setting everything behind back plane to fat
#pragma omp parallel for
  for(int i=0; i<backPlaneInd; i++){
    for(int j=0; j<dim[1]; j++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == innerVal){
	  // set to fat
	  p[0] = ufat;
	}
      }
    }
  }
	
  double boundaryDev = vm["boundary.maxDeviation"].as<double>();

  // resolution of Voronoi segmentation
  double segSize = 0.2;

  int voxSkip = static_cast<int>(floor(segSize/imgRes));
	
  // other side of back plane, do segmentation
#pragma omp parallel for schedule(static,1)
  for(int i=backPlaneInd; i<=dim[0]-voxSkip; i+=voxSkip){
    double coords[3];
    coords[0] = originCoords[0] + imgRes*i;
    for(int j=0; j<=dim[1]-voxSkip; j+=voxSkip){
      coords[1] = originCoords[1] + imgRes*j;
      for(int k=0; k<=dim[2]-voxSkip; k+=voxSkip){

	bool doSeg = false;
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	
	if(p[0] == innerVal){
	  doSeg = true;
	} else {
	  // check other voxels in supervoxel
	  for(int a=0; a<voxSkip; a++){
	    for(int b=0; b<voxSkip; b++){
	      for(int c=0; c<voxSkip; c++){
		p = static_cast<unsigned char*>(breast->GetScalarPointer(i+a,j+b,k+c));
		if(p[0] == innerVal){
		  doSeg = true;
		}
	      }
	    }
	  }
	}
	if(doSeg){
	  // found voxel to segment
	  
	  // find coordinates
	  coords[2] = originCoords[2] + imgRes*k;
	  
	  // nearest fat points
	  vtkSmartPointer<vtkIdList> nearPts =
	    vtkSmartPointer<vtkIdList>::New();
	  
	  findSeed->FindPointsWithinRadius(compSeedRadius, coords, nearPts);

	  // find minimum distance
	  // starting with the first
	  vtkIdType id = nearPts->GetId(0);

	  // compute position in local coordinate system
	  vtkVector3d rvec;
	  vtkVector3d localCoords;
	  for(int m=0; m<3; m++){
	    rvec[m] = coords[m]-fatCompartments[id].pos[m];
	  }
	  for(int m=0; m<3; m++){
	    localCoords[m] = rvec.Dot(fatCompartments[id].axis[m]);
	  }
					
	  // compute distance function
	  double minDist = 0.0;
	  double nextMinDist = 0.0;
	  
	  for(int m=0; m<3; m++){
	    minDist += fatCompartments[id].scale[m]*localCoords[m]*localCoords[m];
	  }
	  minDist = minDist/fatCompartments[id].g;
					
	  vtkIdType closestId = id;
	  vtkIdType nextClosestId = id;

	  int numSeedToCheck = nearPts->GetNumberOfIds();

	  // check other fat points to find minimum distance
	  for(int n=1; n<numSeedToCheck; n++){
	    vtkIdType thisId = nearPts->GetId(n);
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-fatCompartments[thisId].pos[m];
	    }
	    for(int m=0; m<3; m++){
	      localCoords[m] = rvec.Dot(fatCompartments[thisId].axis[m]);
	    }

	    // compute distance
	    double dist = 0.0;
	    for(int m=0; m<3; m++){
	      dist += fatCompartments[thisId].scale[m]*localCoords[m]*localCoords[m];
	    }
	    dist = dist/fatCompartments[thisId].g;
	    
	    if(dist < minDist){
	      nextMinDist = minDist;
	      nextClosestId = closestId;
	      minDist = dist;
	      closestId = thisId;
	    } else if(dist < nextMinDist){
	      nextMinDist = dist;
	      nextClosestId = thisId;
	    }
	  }

	  // now check all gland seeds
	  for(int n=0; n<=numBreastCompartments; n++){
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-glandCompartments[n].pos[m];
	    }
	    for(int m=0; m<3; m++){
	      localCoords[m] = rvec.Dot(glandCompartments[n].axis[m]);
	    }

	    // compute distance
	    double dist = 0.0;
	    for(int m=0; m<3; m++){
	      dist += glandCompartments[n].scale[m]*localCoords[m]*localCoords[m];
	    }
	    dist = dist/glandCompartments[n].g;

	    
	    // glandular compartment so add noise
	    dist += boundaryDev*dist*boundary[n].getNoise(localCoords.Normalized().GetData());
	    
	    if(dist < minDist){
	      nextMinDist = minDist;
	      nextClosestId = closestId;
	      minDist = dist;
	      closestId = numFatSeeds+n;
	    } else if(dist < nextMinDist){
	      nextMinDist = dist;
	      nextClosestId = numFatSeeds+n;
	    }
	  }

	  // set tissue type
	  unsigned char myTissue;
	  
	  if(closestId < numFatSeeds){
	    myTissue = ufat;
	  } else {
	    myTissue = compartmentVal[glandCompartments[closestId-numFatSeeds].compId];
	  }

	  for(int a=0; a<voxSkip; a++){
            for(int b=0; b<voxSkip; b++){
              for(int c=0; c<voxSkip; c++){
                p = static_cast<unsigned char*>(breast->GetScalarPointer(i+a,j+b,k+c));
                if(p[0] == innerVal){
                  p[0] = myTissue;
                }
              }
            }
          }

	  // found two closest ids, set tissue type
	  //if(fabs(nextMinDist-minDist) < ligDistThresh){
	    //// could be ligament
	    //if((closestId < numFatSeeds && nextClosestId >= numFatSeeds) ||
	       //(closestId >= numFatSeeds && nextClosestId < numFatSeeds)){
	      
	      //// gland/fat interface
	      //p[0] = tissue.cooper;
	      //} else {
	      //// assign to fat or compartment
	      //if(closestId < numFatSeeds){
		//p[0] = ufat;
		//} else {
		//// in a compartment
		//p[0] = compartmentVal[glandCompartments[closestId-numFatSeeds].compId];
		//}
	      //}
	    //} else {
	    //// assign to fat or compartment
	    //if(closestId < numFatSeeds){
	      //p[0] = ufat;
	      //} else {
	      //// in a compartment
	      //p[0] = compartmentVal[glandCompartments[closestId-numFatSeeds].compId];
	      //}
	    //}
	}
      }
    }
  }

  // deleting boundary noise
  for(int i=0; i<=numBreastCompartments; i++){
    boundary[i].~perlinNoise();
  }
  ::operator delete(boundary);


  // calculate voxel counts and bounding boxes
  // only updating boundBox for gland compartments
#pragma omp parallel for
  for(int i=0; i<=numBreastCompartments; i++){
    unsigned char val;
    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer());
    for(int c=0; c<dim[2]; c++){
      for(int b=0; b<dim[1]; b++){
	for(int a=0; a<dim[0]; a++){
	  val = *p;
	  p++;
	  if(val == compartmentVal[glandCompartments[i].compId]){
	    if(a < glandCompartments[i].boundBox[0]){
	      glandCompartments[i].boundBox[0] = a;
	    } else {
	      if(a > glandCompartments[i].boundBox[1]){
		glandCompartments[i].boundBox[1] = a;
	      }
	    }
	    if(b < glandCompartments[i].boundBox[2]){
	      glandCompartments[i].boundBox[2] = b;
	    } else {
	      if(b > glandCompartments[i].boundBox[3]){
		glandCompartments[i].boundBox[3] = b;
	      }
	    }
	    if(c < glandCompartments[i].boundBox[4]){
	      glandCompartments[i].boundBox[4] = c;
	    } else {
	      if(c > glandCompartments[i].boundBox[5]){
		glandCompartments[i].boundBox[5] = c;
	      }
	    }
	    glandCompartments[i].voxelCount += 1;
	  }
	}
      }
    }
  }

  // amount of fat and gland and ligament
  double fatVol = 0.0;		// keep track of fat and glandular segmented volume
  double glandVol = 0.0;
  double cooperVol = 0.0;
  double targetGlandVol;
  long long int glandVoxels = 0;
  long long int fatVoxels = 0;
  long long int cooperVoxels = 0;
  //double ligVol = 0.0;

  double targetFatFrac = vm["base.targetFatFrac"].as<double>(); // desired fat fraction of breast
  if(targetFatFrac > 1.0){
    targetFatFrac = 1.0;
  }
  if(targetFatFrac < 0.0){
    targetFatFrac = 0.0;
  }

  // convert gland compartment voxel counts to volume
  //#pragma omp parallel for reduction(+:glandVol,glandVoxels)
  for(int i=0; i<=numBreastCompartments; i++){
    glandCompartments[i].volume = voxelVol*glandCompartments[i].voxelCount;
    glandVol += glandCompartments[i].volume;
    glandVoxels += glandCompartments[i].voxelCount;
  }

  // debug
  //cout << "glandVoxels = " << glandVoxels << "\n";

  // fat and ligament volume need to count voxel by voxel to capture backplane
#pragma omp parallel for reduction(+:fatVoxels)
  for(int c=0; c<dim[2]; c++){
    for(int b=0; b<dim[1]; b++){
      for(int a=0; a<dim[0]; a++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
	if(p[0] == ufat){
	  fatVoxels += 1;
	} else if(p[0] == tissue.cooper){
	  cooperVoxels += 1;
	}
      }
    }
  }
  fatVol = voxelVol*fatVoxels;
  
  //cout << "done.\n";
  //cout << "Initial Voronoi fat fraction = " << fatVol/(glandVol+fatVol) << "\n";


  targetGlandVol = (glandVol+fatVol)*(1-targetFatFrac);

  // save segmented breast
  //vtkSmartPointer<vtkXMLImageDataWriter> writerSeg =
  //	vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //writerSeg->SetFileName(segmentedVoxelFilename.c_str());
  //#if VTK_MAJOR_VERSION <= 5
  //writerSeg->SetInput(breast);
  //#else
  //writerSeg->SetInputData(breast);
  //#endif
  //writerSeg->Write();

  // add back plane voxels assigned to glandular compartment
  int numBackPlaneSkin = 0;
  for(int i=0; i<dim[1]; i++){
    for(int j=0; j<dim[2]; j++){
      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(backPlaneInd,i,j));
      if(p[0] <= compMax && p[0] >= compMin){
        numBackPlaneSkin += 1;
        int ijk[3] =  {backPlaneInd,i,j};
        boundaryList->InsertNextId(breast->ComputePointId(ijk));
      }
    }
  }


  // check boundary voxels and add fat, muscle and near-nipple voxels to delete mask
  vtkIdType nBoundary = boundaryList->GetNumberOfIds();
  vtkIdType remBoundary = nBoundary;

  unsigned char *boundaryDone;
  boundaryDone = new unsigned char [nBoundary];

  for(vtkIdType i=0; i<nBoundary; i++){
    boundaryDone[i] = 0;
  }

  // debug - remove this later
  vtkSmartPointer<vtkIdList> boundaryDel =
    vtkSmartPointer<vtkIdList>::New();
	
  for(vtkIdType i=0; i<nBoundary; i++){
    vtkIdType myId = boundaryList->GetId(i);
    double loc[3];
    double pcoords[3];
    int ijk[3];
    breast->GetPoint(myId, loc);
    breast->ComputeStructuredCoordinates(loc, ijk, pcoords);
		
    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(ijk));
		
    if(p[0] == ufat || p[0] == tissue.muscle || vtkMath::Distance2BetweenPoints(loc, nipplePos) < areolaRad*areolaRad*2){
      // set to remove from list
      //boundaryDel->InsertNextId(myId);
      boundaryDone[i] = 1;
      remBoundary++;
    }
  }
	
  // removals
  //vtkIdType nDel = boundaryDel->GetNumberOfIds();
  
  //for(int i=0; i<nDel; i++){
  //  boundaryList->DeleteId(boundaryDel->GetId(i));
  //}

  // determine glandular compartments to remove

  // find smallest compartment
  unsigned int delNext=0;
  double volNext=glandCompartments[0].volume;

  for(int i=1; i<numBreastCompartments; i++){
    if(glandCompartments[i].volume < volNext){
      delNext = i;
      volNext = glandCompartments[i].volume;
    }
  }

  unsigned int densityClass;
  unsigned int numKeepComp;

  // target fat fraction post compartment removal
  double targetFatFracStep1 = targetFatFrac*0.5;
  double currentFatFrac = fatVol/(fatVol+glandVol);
  
  if(currentFatFrac >= targetFatFracStep1){
    targetFatFracStep1 = currentFatFrac;
  }

  // debug
  // disable gland compartment removal
  //double removeGlandVol = (targetFatFracStep1-currentFatFrac)*(fatVol+glandVol);
  double removeGlandVol = 0.0;

  // debug
  //cout << "removeGlandVol = " << removeGlandVol << "\n";

  if(targetFatFrac < 0.4){
    densityClass = 1;
    // dense breast
  } else if(targetFatFrac < 0.75){
    densityClass = 2;
    // heterogeneous breast
  } else if(targetFatFrac < 0.9){
    densityClass = 3;
    // scattered density
  } else{
    densityClass = 4;
    // fatty breast
  }

  // debug
  //cout << "Density class = " << densityClass << "\n";
  numKeepComp = numBreastCompartments;
  // set numKeepComp based on fat fraction
  //if(targetFatFrac < 0.3){
  //  numKeepComp = numBreastCompartments;
  //  // keep all
  //} else if(targetFatFrac > 0.9){
  //  numKeepComp = 4;
  //  if(numBreastCompartments < 4){
  //    numKeepComp = numBreastCompartments;
  //  }
  //} else {
  //  // linear function of fat fraction
  //  numKeepComp = 2 + static_cast<unsigned int>(round((numBreastCompartments-4)/(0.9-0.3)*(0.9-targetFatFrac)));
  //  if(numBreastCompartments < 4){
  //    numKeepComp = numBreastCompartments;
  //  }
  //}

  // identify compartments to remove
  //unsigned int numDelComp = numBreastCompartments - numKeepComp;
  unsigned int *keepCompList = (unsigned int*)malloc(numBreastCompartments*sizeof(unsigned int));
  unsigned int *delCompList = (unsigned int*)malloc(numBreastCompartments*sizeof(unsigned int));

  unsigned int foundComp = 0;
  unsigned int keepComp = 0;

  //while(foundComp < numDelComp){
  //  // random compartment
  //  unsigned int c = static_cast<unsigned int>(floor(rgen->GetRangeValue(0.0, numBreastCompartments)));
  //  rgen->Next();
  //  if(glandCompartments[c].keep == true){
  //    // found one
  //    delCompList[foundComp] = c;
  //    glandCompartments[c].keep = false;
  //    foundComp++;
  //  }
  //}

  while(removeGlandVol > volNext){
    delCompList[foundComp] = delNext;
    glandCompartments[delNext].keep = false;
    removeGlandVol -= volNext;
    foundComp++;
    //cout << "Removing compartment " << foundComp << " , vol = " << volNext << "\n";
    // update next compartment to delete
    volNext = fatVol+glandVol;
    for(int i=0; i<numBreastCompartments; i++){
      if(glandCompartments[i].volume < volNext && glandCompartments[i].keep == true){
	delNext = i;
	volNext = glandCompartments[i].volume;
      }
    }
  }

  // populate keepCompList
  keepComp = 0;
  for(int i=0; i<numBreastCompartments; i++){
    if(glandCompartments[i].keep == true){
      keepCompList[keepComp] = i;
      keepComp++;
    }
  }
    

  // set removed compartments to fat
#pragma omp parallel for
  for(int i=0; i<foundComp; i++){
    int mc = delCompList[i];
    for(int a=glandCompartments[mc].boundBox[0]; a<=glandCompartments[mc].boundBox[1]; a++){
      for(int b=glandCompartments[mc].boundBox[2]; b<=glandCompartments[mc].boundBox[3]; b++){
	for(int c=glandCompartments[mc].boundBox[4]; c<=glandCompartments[mc].boundBox[5]; c++){
	  unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
	  if(p[0] == compartmentVal[glandCompartments[mc].compId]){
	    p[0] = ufat;
	  }
	}
      }
    }
  } 

  // update fat fraction and voxel counts
  for(int i=0; i<foundComp; i++){
    glandVol -= glandCompartments[delCompList[i]].volume;
    fatVol += glandCompartments[delCompList[i]].volume;
  }

  fatVoxels = static_cast<long long int>(round(fatVol/voxelVol));
  glandVoxels = static_cast<long long int>(round(glandVol/voxelVol));

  //cout << "Compartment removal fat fraction = " << fatVol/(fatVol+glandVol) << "\n"; 
  //cout << "Target fat fraction = " << targetFatFrac << "\n";
  //cout << "fatVoxels = " << fatVoxels << ", glandVoxels = " << glandVoxels << "\n";

  /***********************
	Ducts and TDLUs
  ***********************/

  // timing
  //time(&currentTime);
  //cout << "Time interval 4: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  // create duct network
  //cout << "Creating duct/TDLU network...\n";
	
  // File to store duct locations
  char TDLUlocFilename[128];
  sprintf(TDLUlocFilename, "%s/p_%d.loc", outputDir.c_str(), randSeed); 

  // data structure to store TDLU locations and attributes
  // TDLUs are ovoid
  // attributes are main axis radius, off axis radius, main axis direction (unit vector)
  vtkSmartPointer<vtkPoints> *TDLUloc;
  TDLUloc = new vtkSmartPointer<vtkPoints>[numKeepComp];
	
  vtkSmartPointer<vtkDoubleArray> *TDLUattr;
  TDLUattr = new vtkSmartPointer<vtkDoubleArray>[numKeepComp];
	
  for(int i=0; i<numKeepComp; i++){
    TDLUloc[i] = vtkSmartPointer<vtkPoints>::New();
    TDLUattr[i] = vtkSmartPointer<vtkDoubleArray>::New();
    TDLUattr[i]->SetNumberOfComponents(5);
  }

  // step size for searching for start position is minimum voxel dimension
  double breastSpacing[3];
  breast->GetSpacing(breastSpacing);
  double step = breastSpacing[0];
  for(int i=1; i<2; i++){
    if(breastSpacing[i] < step){
      step = breastSpacing[i];
    }
  }

#pragma omp parallel
  {
#pragma omp for
    for(int i=0; i<keepComp; i++){
      
      //#pragma omp critical
      //{
      //	cout << "Breast compartment " << i << ": " << glandCompartments[keepCompList[i]].voxelCount*imgRes*imgRes*imgRes/1000.0 <<
      //	  " cc (" << glandCompartments[keepCompList[i]].voxelCount << " voxels).\n";
      //}

      // find starting direction
      double sdir[3];

      for(int j=0; j<3; j++){
	sdir[j] = glandCompartments[keepCompList[i]].pos[j] - nipplePos[j];
      }
      // normalize
      double snorm = sqrt(sdir[0]*sdir[0]+sdir[1]*sdir[1]+sdir[2]*sdir[2]);
      for(int j=0; j<3; j++){
	sdir[j] = sdir[j]/snorm;
      }

      // find starting position

      // start from nipple, move towards compartment seed until inside compartment
      double currentPos[3] = {nipplePos[0], nipplePos[1], nipplePos[2]};
      bool inCompartment = false;
      double pcoords[3];
      int indicies[3];
      // distance traveled towards seed point
      double dist = 0.0;

      while(!inCompartment && dist < snorm){
	// find voxel containing current position
	breast->ComputeStructuredCoordinates(currentPos, indicies, pcoords);

	// is the voxel in the compartment?
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(indicies));
	
	if(p[0] == compartmentVal[glandCompartments[keepCompList[i]].compId]){
	  // inside compartment
	  inCompartment = true;
	} else {
	  // outside, step towards seed point
	  for(int j=0; j<3; j++){
	    currentPos[j] = currentPos[j] + step*sdir[j];
	  }
	  dist += step;
	}
      }

      // step a little extra along direction to center of compartment
      for(int j=0; j<3; j++){
	currentPos[j] = currentPos[j] + 5*step*sdir[j];
      }
      dist += 5*step;

      if(!inCompartment){
	cout << "Error, could not find starting duct position\n";
	//	return EXIT_FAILURE;
      }

      // create connector from nipple to start position
      // with cubic spline
		
      double finalDuctRad = 0.2*nippleRad;
		
      vtkSmartPointer<vtkParametricSpline> nipcon = 
	vtkSmartPointer<vtkParametricSpline>::New();
      nipcon->ClosedOff();
      nipcon->SetNumberOfPoints(4);
      nipcon->SetPoint(0, nipplePos[0]+step*nippleNorm[0], nipplePos[1]+step*nippleNorm[1], nipplePos[2]+step*nippleNorm[2]);
      nipcon->SetPoint(1, nipplePos[0], nipplePos[1], nipplePos[2]);
      nipcon->SetPoint(2, currentPos[0], currentPos[1], currentPos[2]);
      nipcon->SetPoint(3, currentPos[0]+step*sdir[0], currentPos[1]+step*sdir[1], currentPos[2]+step*sdir[2]);
		
      nipcon->ParameterizeByLengthOff();
      nipcon->SetLeftConstraint(0);
      nipcon->SetRightConstraint(0);
		
      // step along spline from point 1 to point 2
      double nipconStep = imgRes/6.0/dist;	// make 2*dist/imgRes steps from point 1 to 2
      double initRad = vm["ductTree.initRad"].as<double>();
      double pos[3];
      pos[0] = 1.0/3.0 + nipconStep;	// parametric location of point 1 plus a step
      pos[1] = 0.0;	// these always zero
      pos[2] = 0.0;
		
      while(pos[0] <= 2.0/3.0){
	double myRad = finalDuctRad - (pos[0]-1.0/3.0)*(finalDuctRad-initRad)*3.0;
	double myPos[3];
	double pcoords[3];
	int myPosPix[3];
	nipcon->Evaluate(pos, myPos, NULL);
	// find containing voxel
	breast->ComputeStructuredCoordinates(myPos, myPosPix, pcoords);
	// sphere of radius myRad should be duct - only change values under skin
	int pixRad = (int)ceil(myRad/imgRes);
	int thisPix[3];
	for(int a=-pixRad; a<= pixRad; a++){
	  thisPix[0] = myPosPix[0] + a;
	  for(int b=-pixRad; b<= pixRad; b++){
	    thisPix[1] = myPosPix[1] + b;
	    for(int c=-pixRad; c<= pixRad; c++){
	      thisPix[2] = myPosPix[2] + c;
	      double dist2 = imgRes*imgRes*((double)((thisPix[0]-myPosPix[0])*(thisPix[0]-myPosPix[0])+
						     (thisPix[1]-myPosPix[1])*(thisPix[1]-myPosPix[1])+(thisPix[2]-myPosPix[2])*(thisPix[2]-myPosPix[2])));
	      if(dist2 <= myRad*myRad){
		unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(thisPix));
		unsigned char pval;
#pragma omp atomic read
		pval = *p;
		if(pval != tissue.bg && pval != tissue.skin && pval != tissue.nipple){
#pragma omp atomic write
		  *p = tissue.duct;
		}
	      }
	    }
	  }
	}
	pos[0] += nipconStep;
      }
						  
      // create a seed for duct random number generator
      int seed;
#pragma omp critical (randgen)
      {
	seed = static_cast<int>(round(rgen->GetRangeValue(0.0, 1.0)*2147483648));
	rgen->Next();
      }

      // call duct generation function
      generate_duct(breast, vm, TDLUloc[i], TDLUattr[i], compartmentVal[glandCompartments[keepCompList[i]].compId], 
		    glandCompartments[keepCompList[i]].boundBox, &tissue, currentPos, sdir, seed);
    }
  }

  // write TDLU locations to .loc file
  FILE *TDLUlocFile;
  TDLUlocFile = fopen(TDLUlocFilename, "w");

  for(int i=0; i<keepComp; i++){
    vtkIdType numTDLU = TDLUloc[i]->GetNumberOfPoints();
    double myPoint[3];
    for(int j=0; j<numTDLU; j++){
      TDLUloc[i]->GetPoint(j, myPoint);
      fprintf(TDLUlocFile, "%E,%E,%E\n", myPoint[0], myPoint[1], myPoint[2]);
    }
  }

  fclose(TDLUlocFile);

  // ducts/TDLUs complete
  //cout << "done.\n";

  // save segmented breast with ducts
  //vtkSmartPointer<vtkXMLImageDataWriter> writerSeg2 =
  //	vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //writerSeg2->SetFileName("breast_ductTDLU.vti");
  //#if VTK_MAJOR_VERSION <= 5
  //writerSeg2->SetInput(breast);
  //#else
  //writerSeg2->SetInputData(breast);
  //#endif
  //writerSeg2->Write();

  // timing
  //time(&currentTime);
  //cout << "Time interval 5: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime;

  /**********************
   * 
   *  fat lobules
   * 
   *********************/
	
  // calculate glandular bounding box
  int glandBox[6] = {breastDim[0]+1,-1,breastDim[1]+1,-1,breastDim[2]+1,-1};
  for(int i=0; i<=numBreastCompartments; i++){
    if(glandCompartments[i].keep == true){
      glandBox[0] = (glandCompartments[i].boundBox[0] < glandBox[0]) ? glandCompartments[i].boundBox[0] : glandBox[0];
      glandBox[1] = (glandCompartments[i].boundBox[1] > glandBox[1]) ? glandCompartments[i].boundBox[1] : glandBox[1];
      glandBox[2] = (glandCompartments[i].boundBox[2] < glandBox[2]) ? glandCompartments[i].boundBox[2] : glandBox[2];
      glandBox[3] = (glandCompartments[i].boundBox[3] > glandBox[3]) ? glandCompartments[i].boundBox[3] : glandBox[3];
      glandBox[4] = (glandCompartments[i].boundBox[4] < glandBox[4]) ? glandCompartments[i].boundBox[4] : glandBox[4];
      glandBox[5] = (glandCompartments[i].boundBox[5] > glandBox[5]) ? glandCompartments[i].boundBox[5] : glandBox[5];
    }
  }
	
  // re-label all compartments as gland	
#pragma omp parallel for
  for(int c=glandBox[4]; c<=glandBox[5]; c++){
    for(int b=glandBox[2]; b<=glandBox[3]; b++){
      for(int a=glandBox[0]; a<=glandBox[1]; a++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(a,b,c));
	if(p[0] <= compMax && p[0] >= compMin){
	  // glandular
	  p[0] = ugland;
	}
      }
    }
  }

  currentFatFrac = fatVol/(glandVol+fatVol);
  targetFatFrac = vm["base.targetFatFrac"].as<double>();	// desired fat fraction of breast

  // potentially adjust target fat frac
  if(targetFatFrac<currentFatFrac){
    targetFatFrac = currentFatFrac + (1-currentFatFrac)*0.25;
    cout << "Adjusting target fat fraction before lobulation to " << targetFatFrac << "\n";
  }

  /**********************
   * 
   *  skin boundary fat lobules
   * 
   *********************/

  //double minSkinLobuleAxis = 15.0; // min and max lobule axis length (mm)
  //double maxSkinLobuleAxis = 20.0;
  //double minSkinAxialRatio = 0.4; // min and max lobule axial ratio - was 0.13
  //double maxSkinAxialRatio = 0.5;  // was 0.3
  //double skinLobuleDeflectMax = 0.15; // deflection fraction of pi
  //double skinPerturbMax = 0.15;
  //double skinLigThick = 0.2; // ligament thickness (mm);
  //int maxSkinLobules = 1000;

  double minSkinLobuleAxis; // min and max lobule axis length (mm)
  double maxSkinLobuleAxis;
  double minSkinAxialRatio; // min and max lobule axial ratio - was 0.13
  double maxSkinAxialRatio;  // was 0.3
  double skinLobuleDeflectMax; // deflection fraction of pi
  double skinPerturbMax;
  double targetSkinFatFrac;
  int maxSkinLobules;
  double skinLigThick = vm["lig.thickness"].as<double>(); // ligament thickness (mm);  

  // skin fat lobule Perlin noise coefficients
  double skinLobulePerlinScale = 0.007;
  double skinLobulePerlinLac = 1.8;
  double skinLobulePerlinPers = 0.6;
  int skinLobulePerlinOct = 6;;

  // mega size skin fat lobules
  int numMegaLobules;
  // mega lobule size based on 0.1 of breast volume
  double maxMegaLobuleAxis = pow(0.1*breastVol*0.75/pi/0.75/0.75,1.0/3.0);
  double minMegaLobuleAxis = 0.5*maxMegaLobuleAxis;
  double minMegaAxialRatio = 0.5; 
  double maxMegaAxialRatio = 0.8;  
  double megaLobuleDeflectMax = 0.15; 
  double megaPerturbMax = 0.2;


  //  min/max was 15/20

  switch(densityClass){
  case 1: // dense
    minSkinLobuleAxis = 2.0;
    maxSkinLobuleAxis = 4.0;
    minSkinAxialRatio = 0.5;   
    maxSkinAxialRatio = 0.99;
    skinLobuleDeflectMax = 0.1;
    skinPerturbMax = 0.15;
    maxSkinLobules = 1000;
    numMegaLobules = 1;    
    targetSkinFatFrac = currentFatFrac + (targetFatFrac-currentFatFrac)*0.6;
    break;
  case 2:
    minSkinLobuleAxis = 3.0;
    maxSkinLobuleAxis = 7.0;
    minSkinAxialRatio = 0.5;    
    maxSkinAxialRatio = 0.99;
    skinLobuleDeflectMax = 0.15;
    skinPerturbMax = 0.15;
    maxSkinLobules = 800;
    numMegaLobules = 2;
    targetSkinFatFrac = currentFatFrac + (targetFatFrac-currentFatFrac)*0.6;
    break;
  case 3:
    minSkinLobuleAxis = 8.0;
    maxSkinLobuleAxis = 12.0;
    minSkinAxialRatio = 0.5;   
    maxSkinAxialRatio = 0.99;
    skinLobuleDeflectMax = 0.15;
    skinPerturbMax = 0.15;
    maxSkinLobules = 700;
    numMegaLobules = 4;
    targetSkinFatFrac = currentFatFrac + (targetFatFrac-currentFatFrac)*0.5;
    break;
  case 4:
    minSkinLobuleAxis = 8.0;
    maxSkinLobuleAxis = 12.0;
    minSkinAxialRatio = 0.5;   
    maxSkinAxialRatio = 0.99;
    skinLobuleDeflectMax = 0.15;
    skinPerturbMax = 0.15;
    maxSkinLobules = 500;
    numMegaLobules = 6;
    targetSkinFatFrac = currentFatFrac + (targetFatFrac-currentFatFrac)*0.4;
    break;
  }
    
	
  int numSkinLobules = 0;
  //vtkIdType numGlandBoundary = boundaryList->GetNumberOfIds();
  
  // debug
  //cout << "numGlandBoundary = " << numGlandBoundary << "\n";
  //cout << "target fat frac = " << targetFatFrac << "\n";
  //cout << "Before skin lobules current fat = " << currentFatFrac << ", target for skin = " << targetSkinFatFrac << "\n";
  
  while(numSkinLobules < maxSkinLobules && remBoundary > 0 && currentFatFrac < targetSkinFatFrac){
    
    // debug 
    //cout << "numSkinLobules = " << numSkinLobules << " remBoundary = " << remBoundary << " currentFatFrac = " << currentFatFrac << "\n";

    // pick a random voxel from the boundary list
    int seedVox[3];
    int randVal;
    vtkIdType randId;
    //if(numSkinLobules < numMegaLobules){
    //  randVal = (int)(floor(rgen->GetRangeValue(numGlandBoundary-numBackPlaneSkin, numGlandBoundary)));
    //} else { 
    //  randVal = (int)(floor(rgen->GetRangeValue(0, numGlandBoundary)));
    //}
    //rgen->Next();
    //vtkIdType randId = boundaryList->GetId(randVal);
    
    if(numSkinLobules < numMegaLobules){
      bool foundVox = false;
      while(!foundVox){
	randVal = (int)(floor(rgen->GetRangeValue(nBoundary-numBackPlaneSkin, nBoundary)));
	rgen->Next();
	if(!boundaryDone[randVal]){
	  foundVox = true;
	}
      }
    } else {
      bool foundVox = false;
      while(!foundVox){
	randVal = (int)(floor(rgen->GetRangeValue(0, nBoundary)));
	rgen->Next();
	if(!boundaryDone[randVal]){
          foundVox = true;
	}
      }
    }
    randId = boundaryList->GetId(randVal);

    double loc[3];
    double pcoords[3];
    breast->GetPoint(randId, loc);
    breast->ComputeStructuredCoordinates(loc, seedVox, pcoords);
		
    // found a skin seed voxel, pick point within voxel
    double seed[3];
    for(int i=0; i<3; i++){
      seed[i] = loc[i] + rgen->GetRangeValue(-imgRes/2.0, imgRes/2.0);
      rgen->Next();
    }
		
    // found a seed, pick a size
    double axisLen[3];
    if(numSkinLobules < numMegaLobules){
      axisLen[0] = rgen->GetRangeValue(minMegaLobuleAxis,maxMegaLobuleAxis);
      rgen->Next();
      for(int i=1; i<3; i++){
	axisLen[i] = axisLen[0]*rgen->GetRangeValue(minMegaAxialRatio,maxMegaAxialRatio);
	rgen->Next();
      }
    } else {
      axisLen[0] = rgen->GetRangeValue(minSkinLobuleAxis,maxSkinLobuleAxis);
      rgen->Next();
      for(int i=1; i<3; i++){
	axisLen[i] = axisLen[0]*rgen->GetRangeValue(minSkinAxialRatio,maxSkinAxialRatio);
	rgen->Next();
      }
    }

    // scale factors
    double scaleB = axisLen[1]/axisLen[0];
    double scaleC = axisLen[2]/axisLen[0];
    double A = axisLen[0];
		
    double axis0[3];
    double axis1[3];
    double axis2[3];
		
    // principal direction towards nipple
    for(int j=0; j<3; j++){
      axis0[j] = nipplePos[j] - seed[j];
    }
    vtkMath::Normalize(axis0);
		
    // randomize other axes
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();
		
    // jitter axis direction
    double dtheta,dphi;
    dtheta = rgen->GetRangeValue(0,2*pi);
    rgen->Next();
    if(numSkinLobules < numMegaLobules){
      dphi = rgen->GetRangeValue(0,megaLobuleDeflectMax);
    } else {
      dphi = rgen->GetRangeValue(0,skinLobuleDeflectMax);
    }
    rgen->Next();
    double dr = tan(pi*dphi);

    for(int k=0; k<3; k++){
      axis0[k] = axis0[k] + dr*cos(dtheta)*axis1[k] +
	dr*sin(dtheta)*axis2[k];
    }
		
    vtkMath::Normalize(axis0);
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // debug
    //cout << "Skin lobule dir = " << axis0[0] << " " << axis0[1] << " " << axis0[2] << "\n";

    // convert to vtkVector3d array
    vtkVector3d axis[3];
    for(int i=0; i<3; i++){
      axis[0][i] = axis0[i];
      axis[1][i] = axis1[i];
      axis[2][i] = axis2[i];
    }	

    // Perlin noise for perturbation
    int32_t perturbSeed = (int32_t)(ceil(rgen->GetRangeValue(-1073741824, 1073741824)));
    rgen->Next();
    perlinNoise perturb(vm, perturbSeed, A*skinLobulePerlinScale, skinLobulePerlinLac, skinLobulePerlinPers, skinLobulePerlinOct);
				
    double searchRad;
    if(numSkinLobules < numMegaLobules){
      searchRad = A*(1+megaPerturbMax);
    } else {
      searchRad = A*(1+skinPerturbMax);
    }
    double originalA = A;
		
    // calculate search extent
    int segSpace[6];
    int pixelA = (int)(ceil(searchRad/imgRes));

    segSpace[0] = (seedVox[0] - (int)(pixelA*1.2) > glandBox[0]) ? seedVox[0] - (int)(pixelA*1.2) : glandBox[0];
    segSpace[1] = (seedVox[0] + (int)(pixelA*1.2) < glandBox[1]) ? seedVox[0] + (int)(pixelA*1.2) : glandBox[1];
    segSpace[2] = (seedVox[1] - (int)(pixelA*1.2) > glandBox[2]) ? seedVox[1] - (int)(pixelA*1.2) : glandBox[2];
    segSpace[3] = (seedVox[1] + (int)(pixelA*1.2) < glandBox[3]) ? seedVox[1] + (int)(pixelA*1.2) : glandBox[3];
    segSpace[4] = (seedVox[2] - (int)(pixelA*1.2) > glandBox[4]) ? seedVox[2] - (int)(pixelA*1.2) : glandBox[4];
    segSpace[5] = (seedVox[2] + (int)(pixelA*1.2) < glandBox[5]) ? seedVox[2] + (int)(pixelA*1.2) : glandBox[5];
		
    // iterative over search space, adjusting A as we go
#pragma omp parallel for collapse(3)
    for(int i=segSpace[0]; i<= segSpace[1]; i++){
      for(int j=segSpace[2]; j<= segSpace[3]; j++){
	for(int k=segSpace[4]; k<= segSpace[5]; k++){
	  // continuous standard coordinates
	  double coords[3];
					
	  coords[0] = originCoords[0] + i*imgRes;
	  coords[1] = originCoords[1] + j*imgRes;
	  coords[2] = originCoords[2] + k*imgRes;
					
	  // if duct/TDLU, may need to adjust size
	  unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
					
	  if(*p == tissue.TDLU || *p == tissue.duct){
	    // adjust A
	    
	    // local coordinates
	    vtkVector3d rvec;
	    vtkVector3d lcoords;
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-seed[m];
	    }
	    for(int m=0; m<3; m++){
	      lcoords[m] = rvec.Dot(axis[m]);
	    }
	    
	    // spherical coordinates
	    double r = rvec.Norm();
	    double phi = acos(lcoords[2]/r);
	    double theta = atan2(lcoords[1],lcoords[0]);
						
	    // perturb value based on position on sphere (originalA,phi,theta)
	    double spherePos[3] = {originalA*sin(phi)*cos(theta), originalA*sin(phi)*sin(theta), originalA*cos(phi)};
	    double perturbVal;
	    if(numSkinLobules < numMegaLobules){
	      perturbVal = megaPerturbMax*perturb.getNoise(spherePos);
	    } else {
	      perturbVal = skinPerturbMax*perturb.getNoise(spherePos);
	    }
	    //double bufferVal = 0.5*bufferMax + 0.5*bufferMax*buffer.getNoise(spherePos);
						
	    // in lobule condition is r <= A*(f(theta,phi,scaleB,scaleC)+perturb) if TDLU/duct

	    double f = pow(fabs(sin(phi)*cos(theta)),2.5);
	    f += pow(fabs(sin(phi)*sin(theta))/scaleB,2.5);
	    f += pow(fabs(cos(phi))/scaleC,2.5);
	    f = 1/pow(f,1/2.5);

	    // inside lobule?
	    if(r <= A*(f + perturbVal)){
	      // encroaching, need to adjust A
	      
	      double newA;
	      newA = r/(f + perturbVal);
#pragma omp critical
	      {
		if(A > newA){
		  A = newA;
		}
	      }
	    }
	  }
	}
      }
    }
		
    // have final A, do lobulation
		
    // adjust search radius
    searchRad = searchRad - originalA + A;

    pixelA = (int)(ceil(searchRad/imgRes));

    segSpace[0] = (seedVox[0] - (int)(pixelA*1.2) > glandBox[0]) ? seedVox[0] - (int)(pixelA*1.2) : glandBox[0];
    segSpace[1] = (seedVox[0] + (int)(pixelA*1.2) < glandBox[1]) ? seedVox[0] + (int)(pixelA*1.2) : glandBox[1];
    segSpace[2] = (seedVox[1] - (int)(pixelA*1.2) > glandBox[2]) ? seedVox[1] - (int)(pixelA*1.2) : glandBox[2];
    segSpace[3] = (seedVox[1] + (int)(pixelA*1.2) < glandBox[3]) ? seedVox[1] + (int)(pixelA*1.2) : glandBox[3];
    segSpace[4] = (seedVox[2] - (int)(pixelA*1.2) > glandBox[4]) ? seedVox[2] - (int)(pixelA*1.2) : glandBox[4];
    segSpace[5] = (seedVox[2] + (int)(pixelA*1.2) < glandBox[5]) ? seedVox[2] + (int)(pixelA*1.2) : glandBox[5];

    // iterative over search space, and segment
#pragma omp parallel for collapse(3)
    for(int i=segSpace[0]; i<= segSpace[1]; i++){
      for(int j=segSpace[2]; j<= segSpace[3]; j++){
	for(int k=segSpace[4]; k<= segSpace[5]; k++){
	  
	  // convert glandular tissue
	  unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
					
	  if(*p == ugland){
					
	    // continuous standard coordinates
	    double coords[3];
	    coords[0] = originCoords[0] + i*imgRes;
	    coords[1] = originCoords[1] + j*imgRes;
	    coords[2] = originCoords[2] + k*imgRes;

	    // local coordinates
	    vtkVector3d rvec;
	    vtkVector3d lcoords;
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-seed[m];
	    }
	    for(int m=0; m<3; m++){
	      lcoords[m] = rvec.Dot(axis[m]);
	    }

	    // spherical coordinates
	    double r = rvec.Norm();
	    double phi = acos(lcoords[2]/r);
	    double theta = atan2(lcoords[1],lcoords[0]);

	    // perturb value based on position on sphere (originalA,phi,theta)
	    double spherePos[3] = {originalA*sin(phi)*cos(theta), originalA*sin(phi)*sin(theta), originalA*cos(phi)};
	    double perturbVal;
            if(numSkinLobules < numMegaLobules){
              perturbVal = megaPerturbMax*perturb.getNoise(spherePos);
	    } else {
	      perturbVal = skinPerturbMax*perturb.getNoise(spherePos);
            }
	    
	    double f = pow(fabs(sin(phi)*cos(theta)),2.5);
	    f += pow(fabs(sin(phi)*sin(theta))/scaleB,2.5);
	    f += pow(fabs(cos(phi))/scaleC,2.5);
	    f = 1/pow(f,1/2.5);
						
	    // inside lobule?
	    if(r <= A*(f + perturbVal)-skinLigThick){
	      // gland to fat
	      *p = ufat;
#pragma omp atomic
	      glandVoxels -= 1;
#pragma omp atomic
	      fatVoxels += 1;
	    } else if(r <= A*(f + perturbVal)) {
	      // debug - disabled skin lobule ligaments
	      //	      *p = tissue.cooper;
	      //#pragma omp atomic
	      //cooperVoxels += 1;
	      *p = ufat;
#pragma omp atomic
	      fatVoxels += 1;
#pragma omp atomic
	      glandVoxels -= 1;
	    } 
	  }
	}
      }
    }
		
    // update fatfrac
    currentFatFrac = (double)(fatVoxels)/(double)(fatVoxels+glandVoxels+cooperVoxels);

    glandVol = (double)glandVoxels*pow(imgRes,3.0);
    cooperVol = (double)cooperVoxels*pow(imgRes,3.0);
    fatVol = (double)fatVoxels*pow(imgRes,3.0);
		
    numSkinLobules++;

    // debug
    //cout << "voxel counts g/f/c = " << glandVoxels << " " << fatVoxels << " " << cooperVoxels << "\n";
    
    // update skin boundary
#pragma omp parallel for
    for(int i=0; i<nBoundary; i++){
      if(!boundaryDone[i]){
	// check if should be removed
	vtkIdType myId = boundaryList->GetId(i);
	double loc[3];
	double pcoords[3];
	int ijk[3];
	breast->GetPoint(myId, loc);
	breast->ComputeStructuredCoordinates(loc, ijk, pcoords);
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(ijk));
	
	if(p[0] == ufat || p[0] == tissue.cooper){
	  boundaryDone[i] = 1;
#pragma omp atomic
	  remBoundary--;
	}
      }
    }

    // update skin boundary
    //boundaryDel->Reset();

    /*
    for(int i=0; i<numGlandBoundary; i++){
      vtkIdType myId = boundaryList->GetId(i);
      double loc[3];
      double pcoords[3];
      int ijk[3];
      breast->GetPoint(myId, loc);
      breast->ComputeStructuredCoordinates(loc, ijk, pcoords);
		
      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(ijk));
		
      if(p[0] == ufat || p[0] == tissue.cooper){
	// set to remove from list
	boundaryDel->InsertNextId(myId);
	if(i > numGlandBoundary - numBackPlaneSkin){
	  numBackPlaneSkin--;
	}
      }
    }
	
    // removals
    vtkIdType nDel = boundaryDel->GetNumberOfIds();
    
    for(int i=0; i<nDel; i++){
      boundaryList->DeleteId(boundaryDel->GetId(i));
    }
		
    // update gland boundary	
    numGlandBoundary = boundaryList->GetNumberOfIds();
    */
    // debug
    //if(!(numSkinLobules % 10)){
    //  cout << "skin lobule " << numSkinLobules << ", fat frac = " << currentFatFrac 
    //   << ", target skin ff = " << targetSkinFatFrac << ", numGlandBoundary = " << numGlandBoundary << "\n";
    //}

    // debug
    //cout << "end loop: numSkinLobules = " << numSkinLobules << " remBoundary = " << remBoundary << " currentFatFrac = " << currentFatFrac << "\n";
    //cout << "fatVoxels = " << fatVoxels << " glandVoxels = " << glandVoxels << " cooperVoxels = " << cooperVoxels << "\n";

  }

  // timing
  //time(&currentTime);
  //cout << "Time interval 6: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime;


	
  // debug
  //cout << "Fatfrac after skin lobules = " << currentFatFrac << "\n";

  // save list of points still on skin boundary
  //vtkSmartPointer<vtkPoints> leftPts = 
  //	vtkSmartPointer<vtkPoints>::New();
  
  //for(int i=0; i<boundaryList->GetNumberOfIds(); i++){
  //	vtkIdType myId = boundaryList->GetId(i);
  //	double loc[3];
  //	breast->GetPoint(myId, loc);
  //	leftPts->InsertNextPoint(loc);
  //}
	
  // create a vtkPolyData
  //vtkSmartPointer<vtkPolyData> leftSet =
  //	vtkSmartPointer<vtkPolyData>::New();
  //leftSet->SetPoints(leftPts);

  // add verticies
  //vtkSmartPointer<vtkVertexGlyphFilter> vertexGlyphFilter3 =
  //	vtkSmartPointer<vtkVertexGlyphFilter>::New();
	
  //vertexGlyphFilter3->AddInputData(leftSet);
  //vertexGlyphFilter3->Update();
	
  // output
  //vtkSmartPointer<vtkXMLPolyDataWriter> writerD =
  //	vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  //writerD->SetFileName("leftPoints.vtp");
  //writerD->SetInputConnection(vertexGlyphFilter3->GetOutputPort());
  //#if VTK_MAJOR_VERSION <= 5
  //writerD->SetInput(vertexGlyphFilter3->GetOutput());
  //#else
  //writerD->SetInputConnection(vertexGlyphFilter3->GetOutputPort());
  //#endif
  //writerD->Write();
	
  // save breast with skin lobules
  //vtkSmartPointer<vtkXMLImageDataWriter> writerSegD =
  //	vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //writerSegD->SetFileName("breast_skinlobules.vti");
  //#if VTK_MAJOR_VERSION <= 5
  //writerSegD->SetInput(breast);
  //#else
  //writerSegD->SetInputData(breast);
  //#endif
  //writerSegD->Write();
		
  //return EXIT_SUCCESS;
	
  // finished skin boundary fat lobules

  //cout << "Start internal lobules\n";
	
  /***********************
   * Fat lobules in glandular tissue
   **********************/
	
  //int maxInnerFatLobuleTry = 1600;
  //int numInnerFatLobuleTry = 0;
  //double minInnerLobuleAxis = 3.0;
  //double maxInnerLobuleAxis = 20.0;
  //double innerLobuleDeflectMax = 0.25;
  //double minInnerAxialRatio = 0.1;
  //double maxInnerAxialRatio = 0.5;
  //double innerBufferMax = 0.05;
  //double innerPerturbMax = 0.2;

  int maxInnerFatLobuleTry;
  int numInnerFatLobuleTry = 0;
  double minInnerLobuleAxis;	
  double maxInnerLobuleAxis;
  double innerLobuleDeflectMax;
  double minInnerAxialRatio;
  double maxInnerAxialRatio;
  double innerBufferMax;
  double innerPerturbMax;
  //double innerLigThick = 0.2;

  // inner fat lobule Perlin noise coefficients                                                   
  double innerLobulePerlinScale = 0.007;
  double innerLobulePerlinLac = 1.8;
  double innerLobulePerlinPers = 0.6;
  int innerLobulePerlinOct = 6; 

  // debug
  double minSmallLobuleAxis = 1.0;
  double maxSmallLobuleAxis = 3.0;
  double minSmallAxialRatio = 0.4;
  double maxSmallAxialRatio = 0.8;

  double smallLobuleFraction;


  switch(densityClass){
  case 1: // dense
    maxInnerFatLobuleTry = 1600;
    minInnerLobuleAxis = 2.0;   // was 2/10
    maxInnerLobuleAxis = 5.0;
    innerLobuleDeflectMax = 0.25;
    minInnerAxialRatio = 0.3;  // was 0.2/0.8
    maxInnerAxialRatio = 0.6;
    innerPerturbMax = 0.3;
    smallLobuleFraction = 0.1;
    break;
  case 2: // heterog
    maxInnerFatLobuleTry = 1600;
    minInnerLobuleAxis= 3.0;    // Previous:  3/8 3/15
    maxInnerLobuleAxis= 8.0;
    innerLobuleDeflectMax = 0.25;
    minInnerAxialRatio = 0.2;   // was 0.1/0.5
    maxInnerAxialRatio = 0.6;
    innerPerturbMax = 0.25;
    smallLobuleFraction = 0.15; // Previous: 0.1
    break;
  case 3: // scattered
    maxInnerFatLobuleTry = 1600;
    minInnerLobuleAxis= 3.0;    // was 3/20
    maxInnerLobuleAxis= 12.0;
    innerLobuleDeflectMax = 0.25;
    minInnerAxialRatio = 0.2;     // was 0.1/0.5
    maxInnerAxialRatio = 0.6;
    innerPerturbMax = 0.25;
    smallLobuleFraction = 0.05;
    break;
  case 4: // fatty
    maxInnerFatLobuleTry = 1600;
    minInnerLobuleAxis= 5.0;    // was 5/20
    maxInnerLobuleAxis= 20.0;
    innerLobuleDeflectMax = 0.25;
    minInnerAxialRatio = 0.2;    // was 0.1/0.5
    maxInnerAxialRatio = 0.6;
    innerPerturbMax = 0.2;
    smallLobuleFraction = 0.05;
    break;
  }


  // debug
  //cout << "current ff = " << currentFatFrac << " glandular fat lobule target ff = " << targetFatFrac << "\n";
	
  while(numInnerFatLobuleTry < maxInnerFatLobuleTry && currentFatFrac < targetFatFrac){
   
    // find a seed point in glandular tissue
    bool foundSeed = false;
    // debug
    long long int oldFatVoxels = fatVoxels;
		
    int seedVox[3];
    unsigned char *p;
    while(!foundSeed){
      seedVox[0] = glandBox[0] + (int)(ceil(rgen->GetRangeValue(0, glandBox[1]-glandBox[0])));
      rgen->Next();
      seedVox[1] = glandBox[2] + (int)(ceil(rgen->GetRangeValue(0, glandBox[3]-glandBox[2])));
      rgen->Next();
      seedVox[2] = glandBox[4] + (int)(ceil(rgen->GetRangeValue(0, glandBox[5]-glandBox[4])));
      rgen->Next();
			
      // tissue type
      p = static_cast<unsigned char*>(breast->GetScalarPointer(seedVox));
      if (*p == ugland){
	foundSeed = true;
      }
    }
		
    // found a seed voxel, pick point within voxel
    double seed[3];
    for(int i=0; i<3; i++){
      seed[i] = originCoords[i] + seedVox[i]*imgRes + rgen->GetRangeValue(-imgRes/2.0, imgRes/2.0);
      rgen->Next();
    }
		
    // found a seed, pick a size
    double axisLen[3];
    //axisLen[0] = rgen->GetRangeValue(minInnerLobuleAxis*pow(0.9988,numInnerFatLobuleTry),maxInnerLobuleAxis*pow(0.9988,numInnerFatLobuleTry));
    //axisLen[0] = rgen->GetRangeValue(minInnerLobuleAxis, maxInnerLobuleAxis);
		
    //rgen->Next();
    // debug - base on remaining gland volume to remove
    double remVol = glandVol+cooperVol-targetGlandVol;
    // debug
    //cout << "glandVol = " << glandVol << " cooper = " << cooperVol << " targetGlandVol = " << targetGlandVol << "\n"; 
    if(targetFatFrac-currentFatFrac > smallLobuleFraction){
      axisLen[0] = pow((remVol*0.05+1.0)*0.75/pi/pow((maxInnerAxialRatio+minInnerAxialRatio)/2.0,2.0),1.0/3.0);
      for(int i=1; i<3; i++){
	axisLen[i] = axisLen[0]*rgen->GetRangeValue(minInnerAxialRatio,maxInnerAxialRatio);
	rgen->Next();
      }
    } else {
      axisLen[0] = rgen->GetRangeValue(minSmallLobuleAxis, maxSmallLobuleAxis);
      rgen->Next();
      for(int i=1; i<3; i++){
	axisLen[i] = axisLen[0]*rgen->GetRangeValue(minSmallAxialRatio,maxSmallAxialRatio);
	rgen->Next();
      }
    }
    //cout << "Axis length = " << axisLen[0] << " remVol = " << remVol << "glandVoxels = " << glandVoxels << "\n";


    //for(int i=1; i<3; i++){
    //  axisLen[i] = axisLen[0]*rgen->GetRangeValue(minInnerAxialRatio,maxInnerAxialRatio);
    //  rgen->Next();
    //}

    // scale factors
    double scaleB = axisLen[1]/axisLen[0];
    double scaleC = axisLen[2]/axisLen[0];
    double A = axisLen[0];
		
    double axis0[3];
    double axis1[3];
    double axis2[3];
		
    // principal direction towards nipple
    for(int j=0; j<3; j++){
      axis0[j] = nipplePos[j] - seed[j];
    }
    vtkMath::Normalize(axis0);
		
    // randomize other axes
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // jitter axis direction
    double dtheta,dphi;
    dtheta = rgen->GetRangeValue(0,2*pi);
    rgen->Next();
    dphi = rgen->GetRangeValue(0,innerLobuleDeflectMax);
    rgen->Next();
    double dr = tan(pi*dphi);

    for(int k=0; k<3; k++){
      axis0[k] = axis0[k] + dr*cos(dtheta)*axis1[k] +
	dr*sin(dtheta)*axis2[k];
    }
		
    vtkMath::Normalize(axis0);
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // convert to vtkVector3d array
    vtkVector3d axis[3];
    for(int i=0; i<3; i++){
      axis[0][i] = axis0[i];
      axis[1][i] = axis1[i];
      axis[2][i] = axis2[i];
    }	
		
    // Perlin noise for perturbation and buffer
    int32_t perturbSeed = (int32_t)(ceil(rgen->GetRangeValue(-1073741824, 1073741824)));
    rgen->Next();
    //int32_t bufferSeed = (int32_t)(ceil(rgen->GetRangeValue(-1073741824, 1073741824)));
    //rgen->Next();
    // debug
    //perlinNoise perturb(vm, perturbSeed, "perturb");
    //perlinNoise perturb(vm, perturbSeed, A*0.002, 1.25, 0.75, 6);
    perlinNoise perturb(vm, perturbSeed, A*innerLobulePerlinScale, innerLobulePerlinLac, innerLobulePerlinPers, innerLobulePerlinOct);
    //perlinNoise buffer(vm, bufferSeed, "buffer");
		
    double searchRad = A*(1+innerPerturbMax);
    double originalA = A;
		
    // calculate search extent
    int segSpace[6];
    int pixelA = (int)(ceil(searchRad/imgRes));

    segSpace[0] = (seedVox[0] - (int)(pixelA*1.2) > glandBox[0]) ? seedVox[0] - (int)(pixelA*1.2) : glandBox[0];
    segSpace[1] = (seedVox[0] + (int)(pixelA*1.2) < glandBox[1]) ? seedVox[0] + (int)(pixelA*1.2) : glandBox[1];
    segSpace[2] = (seedVox[1] - (int)(pixelA*1.2) > glandBox[2]) ? seedVox[1] - (int)(pixelA*1.2) : glandBox[2];
    segSpace[3] = (seedVox[1] + (int)(pixelA*1.2) < glandBox[3]) ? seedVox[1] + (int)(pixelA*1.2) : glandBox[3];
    segSpace[4] = (seedVox[2] - (int)(pixelA*1.2) > glandBox[4]) ? seedVox[2] - (int)(pixelA*1.2) : glandBox[4];
    segSpace[5] = (seedVox[2] + (int)(pixelA*1.2) < glandBox[5]) ? seedVox[2] + (int)(pixelA*1.2) : glandBox[5];
		
    // iterative over search space, adjusting A as we go
    /*
      #pragma omp parallel for collapse(3)
      for(int i=segSpace[0]; i<= segSpace[1]; i++){
      for(int j=segSpace[2]; j<= segSpace[3]; j++){
      for(int k=segSpace[4]; k<= segSpace[5]; k++){
      // continuous standard coordinates
      double coords[3];
					
      coords[0] = originCoords[0] + i*imgRes;
      coords[1] = originCoords[1] + j*imgRes;
      coords[2] = originCoords[2] + k*imgRes;
					
      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
      
      // shrink if hit buffer
      if(*p == glandBuffer){
      // adjust A if currently inside lobule
      
      // local coordinates
      vtkVector3d rvec;
      vtkVector3d lcoords;
      for(int m=0; m<3; m++){
      rvec[m] = coords[m]-seed[m];
      }
      for(int m=0; m<3; m++){
      lcoords[m] = rvec.Dot(axis[m]);
      }
      
      // spherical coordinates
      double r = rvec.Norm();
      double phi = acos(lcoords[2]/r);
      double theta = atan2(lcoords[1],lcoords[0]);
      
      // perturb value based on position on sphere (originalA,phi,theta)
      double spherePos[3] = {originalA*sin(phi)*cos(theta), originalA*sin(phi)*sin(theta), originalA*cos(phi)};
      double perturbVal = innerPerturbMax*perturb.getNoise(spherePos);
      double bufferVal = 0.5*innerBufferMax + 0.5*innerBufferMax*buffer.getNoise(spherePos);
						
      // in lobule condition is r <= A*(f(theta,phi,scaleB,scaleC)+perturb+buffer) if TDLU/duct
      // or r <= A*(f(theta,phi,scaleB,scaleC)+perturb) if buffer

      double f = pow(fabs(sin(phi)*cos(theta)),2.5);
      f += pow(fabs(sin(phi)*sin(theta))/scaleB,2.5);
      f += pow(fabs(cos(phi))/scaleC,2.5);
      f = 1/pow(f,1/2.5);

      // inside lobule?
      if(r < A*(f + perturbVal)){
      // encroaching, need to adjust A
							
      double newA;
				
      newA = r/(f + perturbVal);
      #pragma omp critical
      {
      if(A > newA){
      A = newA;
      }
      }
      }
      }
      }
      }
      }
      
    */
		
    // have final A, do lobulation
		
    // adjust search radius
    //searchRad = searchRad - originalA + A;

    pixelA = (int)(ceil(searchRad/imgRes));

    segSpace[0] = (seedVox[0] - (int)(pixelA*1.2) > glandBox[0]) ? seedVox[0] - (int)(pixelA*1.2) : glandBox[0];
    segSpace[1] = (seedVox[0] + (int)(pixelA*1.2) < glandBox[1]) ? seedVox[0] + (int)(pixelA*1.2) : glandBox[1];
    segSpace[2] = (seedVox[1] - (int)(pixelA*1.2) > glandBox[2]) ? seedVox[1] - (int)(pixelA*1.2) : glandBox[2];
    segSpace[3] = (seedVox[1] + (int)(pixelA*1.2) < glandBox[3]) ? seedVox[1] + (int)(pixelA*1.2) : glandBox[3];
    segSpace[4] = (seedVox[2] - (int)(pixelA*1.2) > glandBox[4]) ? seedVox[2] - (int)(pixelA*1.2) : glandBox[4];
    segSpace[5] = (seedVox[2] + (int)(pixelA*1.2) < glandBox[5]) ? seedVox[2] + (int)(pixelA*1.2) : glandBox[5];

    // iterative over search space, and segment
#pragma omp parallel for collapse(3)
    for(int i=segSpace[0]; i<= segSpace[1]; i++){
      for(int j=segSpace[2]; j<= segSpace[3]; j++){
	for(int k=segSpace[4]; k<= segSpace[5]; k++){
	  
	  // convert glandular tissue
	  unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
					
	  // debug - add overwrite of TDLU
	  if(*p == ugland || *p == tissue.TDLU){
	    
	    // continuous standard coordinates
	    double coords[3];
	    coords[0] = originCoords[0] + i*imgRes;
	    coords[1] = originCoords[1] + j*imgRes;
	    coords[2] = originCoords[2] + k*imgRes;

	    // local coordinates
	    vtkVector3d rvec;
	    vtkVector3d lcoords;
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-seed[m];
	    }
	    for(int m=0; m<3; m++){
	      lcoords[m] = rvec.Dot(axis[m]);
	    }

	    // spherical coordinates
	    double r = rvec.Norm();
	    double phi = acos(lcoords[2]/r);
	    double theta = atan2(lcoords[1],lcoords[0]);

	    // perturb value based on position on sphere (originalA,phi,theta)
	    double spherePos[3] = {originalA*sin(phi)*cos(theta), originalA*sin(phi)*sin(theta), originalA*cos(phi)};
	    double perturbVal = innerPerturbMax*perturb.getNoise(spherePos);
	    //double bufferVal = 0.5*innerBufferMax + 0.5*innerBufferMax*buffer.getNoise(spherePos);
						
	    // in lobule condition is r <= A*(f(theta,phi,scaleB,scaleC)+perturb+buffer) if TDLU/duct
	    // or r <= A*(f(theta,phi,scaleB,scaleC)+perturb) if buffer

	    double f = pow(fabs(sin(phi)*cos(theta)),2.5);
	    f += pow(fabs(sin(phi)*sin(theta))/scaleB,2.5);
	    f += pow(fabs(cos(phi))/scaleC,2.5);
	    f = 1/pow(f,1/2.5);
						
	    // inside lobule?
	    if(r <= A*(f + perturbVal)){
	      // gland to fat
	      *p = ufat;
#pragma omp atomic
	      glandVoxels -= 1;
#pragma omp atomic
	      fatVoxels += 1;
	    }
	  }
	}
      }
    }
    // update fatfrac
    currentFatFrac = (double)(fatVoxels)/(double)(fatVoxels+glandVoxels);
    //if(!(numInnerFatLobuleTry % 10)){
    //  cout << "Finished inner lobule " << numInnerFatLobuleTry+1 << ", fat fraction = " << currentFatFrac << "\n";
    //}
    glandVol = glandVoxels*pow(imgRes,3.0);
    numInnerFatLobuleTry += 1;
  }

  // timing
  //time(&currentTime);
  //cout << "Time interval 7: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime;



	
  // debug
	
  // save segmented breast with ducts
  //vtkSmartPointer<vtkXMLImageDataWriter> writerDebug2 =
  //	vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //writerDebug2->SetFileName("breast_lobulated.vti");
  //#if VTK_MAJOR_VERSION <= 5
  //writerDebug2->SetInput(breast);
  //#else
  //writerDebug2->SetInputData(breast);
  //#endif
  //writerDebug2->Write();
	
  /*******************
   * Cooper's Ligaments
   ******************/

  double ligThick = vm["lig.thickness"].as<double>();
	 
  //const unsigned char ligFat = 99;
  double maxLigAxis = vm["lig.maxAxis"].as<double>();
  double minLigAxis = vm["lig.minAxis"].as<double>();
  double maxLigAxisRatio = vm["lig.maxAxialRatio"].as<double>();
  double minLigAxisRatio = vm["lig.minAxialRatio"].as<double>();
  double ligPerturbMax = vm["lig.maxPerturb"].as<double>();
  double ligDeflectMax = vm["lig.maxDeflect"].as<double>();
  //unsigned int numLigs = 5*(unsigned int)ceil(fatVoxels*voxelVol/(4.0/3.0*pi*pow((maxFatLigAxis+minFatLigAxis)/2.0,3.0))); // was 7

  double ligamentedFrac = 0.0;
  double targetLigFrac = vm["lig.targetFrac"].as<double>();
  long long int ligedVoxels = 0;

  double ligPerlinScale = vm["lig.scale"].as<double>();
  double ligPerlinLac = vm["lig.lacunarity"].as<double>();
  double ligPerlinPers = vm["lig.persistence"].as<double>();
  int ligPerlinOct = vm["lig.numOctaves"].as<int>();
	
  int maxFatLigTry = vm["lig.maxTry"].as<int>();
  int flnum = 0;
  int fltry = 0;

  while(ligamentedFrac < targetLigFrac && fltry < maxFatLigTry){  
		 
    // find a voxel in unligamented region
    bool foundSeed = false;
		
    int seedVox[3];
    unsigned char *p;
    while(!foundSeed){
      seedVox[0] = breastExtent[0] + (int)(ceil(rgen->GetRangeValue(0, breastExtent[1]-breastExtent[0])));
      rgen->Next();
      seedVox[1] = breastExtent[2] + (int)(ceil(rgen->GetRangeValue(0, breastExtent[3]-breastExtent[2])));
      rgen->Next();
      seedVox[2] = breastExtent[4] + (int)(ceil(rgen->GetRangeValue(0, breastExtent[5]-breastExtent[4])));
      rgen->Next();
			
      fltry += 1;
			
      // tissue type
      p = static_cast<unsigned char*>(breast->GetScalarPointer(seedVox));
      if (*p == ufat || *p == ugland){
	foundSeed = true;
      }
    }
		
    // found seed
    flnum += 1;
		
    // found a seed voxel, pick point within voxel
    double seed[3];
    for(int i=0; i<3; i++){
      seed[i] = originCoords[i] + seedVox[i]*imgRes + rgen->GetRangeValue(-imgRes/2.0, imgRes/2.0);
      rgen->Next();
    }
		
    // found a seed, pick a size
    double axisLen[3];
    axisLen[0] = rgen->GetRangeValue(minLigAxis,maxLigAxis);
    rgen->Next();
    for(int i=1; i<3; i++){
      axisLen[i] = axisLen[0]*rgen->GetRangeValue(minLigAxisRatio,maxLigAxisRatio);
      rgen->Next();
    }

    // scale factors
    double scaleB = axisLen[1]/axisLen[0];
    double scaleC = axisLen[2]/axisLen[0];
    double A = axisLen[0];
		
    double axis0[3];
    double axis1[3];
    double axis2[3];
		
    // principal direction towards nearest surface		
    // find closest surface point
    //double nearPos[3];		
    //innerPoly->GetPoint(locator->FindClosestPoint(seed), nearPos);
		
    //for(int j=0; j<3; j++){
    //  axis0[j] = nearPos[j] - seed[j];
    //}
    
    // principle direction towards nipple
    for(int j=0; j<3; j++){
      axis0[j] = nipplePos[j] - seed[j];
    }


    vtkMath::Normalize(axis0);
		
    // randomize other axes
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // jitter axis direction
    double dtheta,dphi;
    dtheta = rgen->GetRangeValue(0,2*pi);
    rgen->Next();
    dphi = rgen->GetRangeValue(0,ligDeflectMax);
    rgen->Next();
    double dr = tan(pi*dphi);

    for(int k=0; k<3; k++){
      axis0[k] = axis0[k] + dr*cos(dtheta)*axis1[k] +
	dr*sin(dtheta)*axis2[k];
    }
		
    vtkMath::Normalize(axis0);
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // convert to vtkVector3d array
    vtkVector3d axis[3];
    for(int i=0; i<3; i++){
      axis[0][i] = axis0[i];
      axis[1][i] = axis1[i];
      axis[2][i] = axis2[i];
    }	

    // TODO: jitter axis directions
    
    // Perlin noise for perturbation and buffer
    int32_t perturbSeed = (int32_t)(ceil(rgen->GetRangeValue(-1073741824, 1073741824)));
    rgen->Next();
    perlinNoise perturb(vm, perturbSeed, A*ligPerlinScale, ligPerlinLac, ligPerlinPers, ligPerlinOct);
				
    double searchRad = A*(1+ligPerturbMax)+ligThick;
		
    // calculate search extent
    int segSpace[6];
    int pixelA = (int)(ceil(searchRad/imgRes));

    segSpace[0] = (seedVox[0] - (int)(pixelA*1.2) > breastExtent[0]) ? seedVox[0] - (int)(pixelA*1.2) : breastExtent[0];
    segSpace[1] = (seedVox[0] + (int)(pixelA*1.2) < breastExtent[1]) ? seedVox[0] + (int)(pixelA*1.2) : breastExtent[1];
    segSpace[2] = (seedVox[1] - (int)(pixelA*1.2) > breastExtent[2]) ? seedVox[1] - (int)(pixelA*1.2) : breastExtent[2];
    segSpace[3] = (seedVox[1] + (int)(pixelA*1.2) < breastExtent[3]) ? seedVox[1] + (int)(pixelA*1.2) : breastExtent[3];
    segSpace[4] = (seedVox[2] - (int)(pixelA*1.2) > breastExtent[4]) ? seedVox[2] - (int)(pixelA*1.2) : breastExtent[4];
    segSpace[5] = (seedVox[2] + (int)(pixelA*1.2) < breastExtent[5]) ? seedVox[2] + (int)(pixelA*1.2) : breastExtent[5];
		
    // iterative over search space, and segment
#pragma omp parallel for collapse(3)
    for(int i=segSpace[0]; i<= segSpace[1]; i++){
      for(int j=segSpace[2]; j<= segSpace[3]; j++){
	for(int k=segSpace[4]; k<= segSpace[5]; k++){
	  
	  // convert glandular tissue
	  unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	  
	  if(*p == ufat || *p == ugland){
					
	    // continuous standard coordinates
	    double coords[3];
	    coords[0] = originCoords[0] + i*imgRes;
	    coords[1] = originCoords[1] + j*imgRes;
	    coords[2] = originCoords[2] + k*imgRes;

	    // local coordinates
	    vtkVector3d rvec;
	    vtkVector3d lcoords;
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-seed[m];
	    }
	    for(int m=0; m<3; m++){
	      lcoords[m] = rvec.Dot(axis[m]);
	    }

	    // spherical coordinates
	    double r = rvec.Norm();
	    double phi = acos(lcoords[2]/r);
	    double theta = atan2(lcoords[1],lcoords[0]);

	    double spherePos[3] = {A*sin(phi)*cos(theta), A*sin(phi)*sin(theta), A*cos(phi)};
	    double perturbVal = ligPerturbMax*perturb.getNoise(spherePos);

	    double f = pow(fabs(sin(phi)*cos(theta)),2.7);
	    f += pow(fabs(sin(phi)*sin(theta))/scaleB,2.7);
	    f += pow(fabs(cos(phi))/scaleC,2.7);
	    f = 1.0/pow(f,1/2.7);
	    
	    // inside ligament lobule?
	    if(r <= A*(f + perturbVal)-ligThick){
	      // interior of ligament volume
	      if(*p == ufat){
		*p = tissue.fat;
	      } else {
		*p = tissue.gland;
	      }
#pragma omp atomic
	      ligedVoxels += 1;
	    } else if(r <= A*(f + perturbVal)) {
	      *p = tissue.cooper;
#pragma omp atomic
	      ligedVoxels += 1;
	    } 
	  }
	}
      }
    }
    // update ligamented frac
    ligamentedFrac = static_cast<double>(ligedVoxels)/static_cast<double>(fatVoxels+glandVoxels);
    //cout << "Lig " << fltry << ", Ligamented fraction = " << ligamentedFrac << "\n";
  }

  // timing
  //time(&currentTime);
  //cout << "Time interval 8: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime;


	
  //cout << "Final ligament fraction = " << ligamentedFrac << "\n";
	
  /*******************
   * Cooper's Ligaments in Glandular Tissue
   ******************/
  /*	 
  const unsigned char ligGland = 160;
  double maxGlandLigAxis = 10.0;
  double minGlandLigAxis = 6.0;
  double maxGlandLigAxisRatio = 0.6;
  double minGlandLigAxisRatio = 0.3;
  double glandLigPerturbMax = 0.3;
  double innerLigDeflectMax = 0.15;
  unsigned int numGlandLigs = 5.0*(unsigned int)ceil(glandVoxels*voxelVol/(4.0/3.0*pi*pow((maxGlandLigAxis+minGlandLigAxis)/2.0,3.0))); // was 3
  int maxGlandLigTry = 10000;
  int glnum = 0;
  int gltry = 0;
  while(glnum < numGlandLigs && gltry < maxGlandLigTry){  
    //for(int lnum=0; lnum<numGlandLigs; lnum++){
    // find a voxel in unligamented fat
    bool foundSeed = false;
    
    int seedVox[3];
    unsigned char *p;
    while(!foundSeed){
      seedVox[0] = glandBox[0] + (int)(ceil(rgen->GetRangeValue(0, glandBox[1]-glandBox[0])));
      rgen->Next();
      seedVox[1] = glandBox[2] + (int)(ceil(rgen->GetRangeValue(0, glandBox[3]-glandBox[2])));
      rgen->Next();
      seedVox[2] = glandBox[4] + (int)(ceil(rgen->GetRangeValue(0, glandBox[5]-glandBox[4])));
      rgen->Next();
      
      gltry += 1;
			
      // tissue type
      p = static_cast<unsigned char*>(breast->GetScalarPointer(seedVox));
      if (*p == tissue.gland || *p == glandBuffer){
	foundSeed = true;
      }
    }
		
    // found seed voxel
    glnum += 1;
		
    // found a seed voxel, pick point within voxel
    double seed[3];
    for(int i=0; i<3; i++){
      seed[i] = originCoords[i] + seedVox[i]*imgRes + rgen->GetRangeValue(-imgRes/2.0, imgRes/2.0);
      rgen->Next();
    }
    
    // found a seed, pick a size
    double axisLen[3];
    axisLen[0] = rgen->GetRangeValue(minGlandLigAxis,maxGlandLigAxis);
    rgen->Next();
    for(int i=1; i<3; i++){
      axisLen[i] = axisLen[0]*rgen->GetRangeValue(minGlandLigAxisRatio,maxGlandLigAxisRatio);
      rgen->Next();
    }

    // scale factors
    double scaleB = axisLen[1]/axisLen[0];
    double scaleC = axisLen[2]/axisLen[0];
    double A = axisLen[0];
		
    double axis0[3];
    double axis1[3];
    double axis2[3];
		
    // principal direction towards nipple
    for(int j=0; j<3; j++){
      axis0[j] = nipplePos[j] - seed[j];
    }
    vtkMath::Normalize(axis0);
		
    // randomize other axes
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // jitter axis direction
    double dtheta,dphi;
    dtheta = rgen->GetRangeValue(0,2*pi);
    rgen->Next();
    dphi = rgen->GetRangeValue(0,innerLigDeflectMax);
    rgen->Next();
    double dr = tan(pi*dphi);

    for(int k=0; k<3; k++){
      axis0[k] = axis0[k] + dr*cos(dtheta)*axis1[k] +
	dr*sin(dtheta)*axis2[k];
    }
		
    vtkMath::Normalize(axis0);
    vtkMath::Perpendiculars(axis0, axis1, axis2, rgen->GetRangeValue(0,2*pi));
    rgen->Next();

    // convert to vtkVector3d array
    vtkVector3d axis[3];
    for(int i=0; i<3; i++){
      axis[0][i] = axis0[i];
      axis[1][i] = axis1[i];
      axis[2][i] = axis2[i];
    }	
		
    // Perlin noise for perturbation
    int32_t perturbSeed = (int32_t)(ceil(rgen->GetRangeValue(-1073741824, 1073741824)));
    rgen->Next();
    perlinNoise perturb(vm, perturbSeed, "perturb");
		
    double searchRad = A*(1+glandLigPerturbMax)+ligThick;
		
    // calculate search extent
    int segSpace[6];
    int pixelA = (int)(ceil(searchRad/imgRes));

    segSpace[0] = (seedVox[0] - (int)(pixelA*1.2) > glandBox[0]) ? seedVox[0] - (int)(pixelA*1.2) : glandBox[0];
    segSpace[1] = (seedVox[0] + (int)(pixelA*1.2) < glandBox[1]) ? seedVox[0] + (int)(pixelA*1.2) : glandBox[1];
    segSpace[2] = (seedVox[1] - (int)(pixelA*1.2) > glandBox[2]) ? seedVox[1] - (int)(pixelA*1.2) : glandBox[2];
    segSpace[3] = (seedVox[1] + (int)(pixelA*1.2) < glandBox[3]) ? seedVox[1] + (int)(pixelA*1.2) : glandBox[3];
    segSpace[4] = (seedVox[2] - (int)(pixelA*1.2) > glandBox[4]) ? seedVox[2] - (int)(pixelA*1.2) : glandBox[4];
    segSpace[5] = (seedVox[2] + (int)(pixelA*1.2) < glandBox[5]) ? seedVox[2] + (int)(pixelA*1.2) : glandBox[5];
		
    // iterative over search space, and segment
#pragma omp parallel for collapse(3)
    for(int i=segSpace[0]; i<= segSpace[1]; i++){
      for(int j=segSpace[2]; j<= segSpace[3]; j++){
	for(int k=segSpace[4]; k<= segSpace[5]; k++){
	  
	  unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
					
	  if(*p == tissue.gland || *p == glandBuffer){
	    
	    // continuous standard coordinates
	    double coords[3];
	    coords[0] = originCoords[0] + i*imgRes;
	    coords[1] = originCoords[1] + j*imgRes;
	    coords[2] = originCoords[2] + k*imgRes;

	    // local coordinates
	    vtkVector3d rvec;
	    vtkVector3d lcoords;
	    for(int m=0; m<3; m++){
	      rvec[m] = coords[m]-seed[m];
	    }
	    for(int m=0; m<3; m++){
	      lcoords[m] = rvec.Dot(axis[m]);
	    }

	    // spherical coordinates
	    double r = rvec.Norm();
	    double phi = acos(lcoords[2]/r);
	    double theta = atan2(lcoords[1],lcoords[0]);
	    
	    double spherePos[3] = {A*sin(phi)*cos(theta), A*sin(phi)*sin(theta), A*cos(phi)};
	    double perturbVal = fatLigPerturbMax*perturb.getNoise(spherePos);
						
	    double f = pow(fabs(sin(phi)*cos(theta)),2.7);
	    f += pow(fabs(sin(phi)*sin(theta))/scaleB,2.7);
	    f += pow(fabs(cos(phi))/scaleC,2.7);
	    //f = 1/sqrt(f);
	    f = 1.0/pow(f,1/2.7);
						
	    // inside ligament lobule?
	    if(r <= A*(f + perturbVal)-ligThick){
	      *p = ligGland;
	    } else if(r <= A*(f + perturbVal)){
	      *p = tissue.cooper;
	    } 
	  }
	}
      }
    }
  }
	
  cout << "Created " << glnum << " of " << numGlandLigs <<" gland ligaments\n";
  */
	
  // convert remaining ufat and ugland 
	
#pragma omp parallel for schedule(static,1)
  for(int k=0; k<dim[2]; k++){
    unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(0,0,k));
    for(int j=0; j<dim[1]; j++){
      for(int i=0; i<dim[0]; i++){
	if(*p == ugland){
	  *p = tissue.gland;
	} else if(*p == ufat){
	  *p = tissue.fat;
	}
	p++;
      }
    }
  }
		
  // calculate gland and fat bounding boxes
  int fatVoxBound[6] = {breastDim[0]+1,-1,breastDim[1]+1,-1,breastDim[2],-1};
  int glandVoxBound[6] = {breastDim[0]+1,-1,breastDim[1]+1,-1,breastDim[2],-1};
	
  for(int i=0; i<dim[0]; i++){
    for(int j=0; j<dim[1]; j++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.fat){
	  fatVoxBound[0] = i;
	  goto foundf0;
	}
      }
    }
  }

 foundf0:
  for(int i=dim[0]-1; i>=0; i--){
    for(int j=0; j<dim[1]; j++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.fat){
	  fatVoxBound[1] = i;
	  goto foundf1;
	}
      }
    }
  }
	
 foundf1:
  for(int j=0; j<dim[1]; j++){
    for(int i=0; i<dim[0]; i++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.fat){
	  fatVoxBound[2] = j;
	  goto foundf2;
	}
      }
    }
  }
	
 foundf2:
  for(int j=dim[1]-1; j>=0; j--){
    for(int i=0; i<dim[0]; i++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.fat){
	  fatVoxBound[3] = j;
	  goto foundf3;
	}
      }
    }
  }		
					
 foundf3:
  for(int k=0; k<dim[2]; k++){
    for(int i=0; i<dim[0]; i++){
      for(int j=0; j<dim[1]; j++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.fat){
	  fatVoxBound[4] = k;
	  goto foundf4;
	}
      }
    }
  }
	
 foundf4:
  for(int k=dim[2]-1; k>=0; k--){
    for(int i=0; i<dim[0]; i++){
      for(int j=0; j<dim[1]; j++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.fat){
	  fatVoxBound[5] = k;
	  goto foundf5;
	}
      }
    }
  }
	
 foundf5:
  for(int i=0; i<dim[0]; i++){
    for(int j=0; j<dim[1]; j++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.duct || p[0] == tissue.TDLU || p[0] == tissue.gland){
	  glandVoxBound[0] = i;
	  goto foundg0;
	}
      }
    }
  }

 foundg0:
  for(int i=dim[0]-1; i>=0; i--){
    for(int j=0; j<dim[1]; j++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.duct || p[0] == tissue.TDLU || p[0] == tissue.gland){
	  glandVoxBound[1] = i;
	  goto foundg1;
	}
      }
    }
  }
	
 foundg1:
  for(int j=0; j<dim[1]; j++){
    for(int i=0; i<dim[0]; i++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.duct || p[0] == tissue.TDLU || p[0] == tissue.gland){
	  glandVoxBound[2] = j;
	  goto foundg2;
	}
      }
    }
  }
	
 foundg2:
  for(int j=dim[1]-1; j>=0; j--){
    for(int i=0; i<dim[0]; i++){
      for(int k=0; k<dim[2]; k++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.duct || p[0] == tissue.TDLU || p[0] == tissue.gland){
	  glandVoxBound[3] = j;
	  goto foundg3;
	}
      }
    }
  }		
					
 foundg3:
  for(int k=0; k<dim[2]; k++){
    for(int i=0; i<dim[0]; i++){
      for(int j=0; j<dim[1]; j++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.duct || p[0] == tissue.TDLU || p[0] == tissue.gland){
	  glandVoxBound[4] = k;
	  goto foundg4;
	}
      }
    }
  }
	
 foundg4:
  for(int k=dim[2]-1; k>=0; k--){
    for(int i=0; i<dim[0]; i++){
      for(int j=0; j<dim[1]; j++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.duct || p[0] == tissue.TDLU || p[0] == tissue.gland){
	  glandVoxBound[5] = k;
	  goto foundg5;
	}
      }
    }
  }
		
  // save voronoi breast with ducts
 foundg5:
  //vtkSmartPointer<vtkXMLImageDataWriter> writerVoron =
  //vtkSmartPointer<vtkXMLImageDataWriter>::New();

  //writerVoron->SetFileName("breastSegmented.vti");
  //#if VTK_MAJOR_VERSION <= 5
  //writerVoron->SetInput(breast);
  //#else
  //writerVoron->SetInputData(breast);
  //#endif
  //writerVoron->Write();

  /********************
   * Vascular network
   *******************/

  // timing
  //time(&currentTime);
  //cout << "Time interval 9: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime;



  // debug
  //time(&currentTime);
  //cout << "Time interval 18: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  // create arterial network
  //cout << "Creating vascular network...\n";

  // arterial sources - internal thoracic (2 -main), thoracocranial, lateral thoracic
  // veinous drainage - auxillary vein (2 - main), subclavian, intercostal, internal thoracic veins
	
  // find breast extent
  int internalExtentVox[6];

  internalExtentVox[0] = (fatVoxBound[0]<glandVoxBound[0]) ? fatVoxBound[0] : glandVoxBound[0];
  internalExtentVox[1] = (fatVoxBound[1]>glandVoxBound[1]) ? fatVoxBound[1] : glandVoxBound[1];
  internalExtentVox[2] = (fatVoxBound[2]<glandVoxBound[2]) ? fatVoxBound[2] : glandVoxBound[2];
  internalExtentVox[3] = (fatVoxBound[3]>glandVoxBound[3]) ? fatVoxBound[3] : glandVoxBound[3];
  internalExtentVox[4] = (fatVoxBound[4]<glandVoxBound[4]) ? fatVoxBound[4] : glandVoxBound[4];
  internalExtentVox[5] = (fatVoxBound[5]>glandVoxBound[5]) ? fatVoxBound[5] : glandVoxBound[5];
	
  int backMinInd = internalExtentVox[0];
	
  // move into breast by 1 mm
  backMinInd += (int)(ceil(1.0/imgRes));
	
  // create mask for back plane
  vtkSmartPointer<vtkImageData> backPlane = 
    vtkSmartPointer<vtkImageData>::New();
		
  backPlane->SetSpacing(breast->GetSpacing());
  
  backPlane->SetExtent(backMinInd, backMinInd, 0, dim[1]-1, 0, dim[2]-1);
  backPlane->SetOrigin(breast->GetOrigin());

  // allocate unsigned char
#if VTK_MAJOR_VERSION <= 5
  backPlane->SetNumberOfScalarComponents(1);
  backPlane->SetScalarTypeToUnsignedChar();
  backPlane->AllocateScalars();
#else
  backPlane->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
		
  double backMass[2] = {0.0, 0.0};
  long long int voxelCount = 0;
	
  // create mask and calculate 2D center of mass
  for(int j=0; j<dim[1]; j++){
    for(int k=0; k<dim[2]; k++){
      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(backMinInd,j,k));
      unsigned char* q = static_cast<unsigned char*>(backPlane->GetScalarPointer(backMinInd,j,k));
      if(p[0] != tissue.bg && p[0] != tissue.skin){
	q[0] = 1;
	voxelCount += 1;
	backMass[0] = backMass[0] + (j - backMass[0])/voxelCount;
	backMass[1] = backMass[1] + (k - backMass[1])/voxelCount;
      } else {
	q[0] = 0;
      }
    }
  }
	
  int backMassVox[3];
  backMassVox[0] = backMinInd;
  backMassVox[1] = (int)(round(backMass[0]));
  backMassVox[2] = (int)(round(backMass[1]));
	
  // calculate continuous position of mass center indicies [backMinInd, backMassVox[0], backMassVox[1]]
  double backCenter[3];
  breast->GetPoint(breast->ComputePointId(backMassVox),backCenter);
	
  double arteryAngle[4];
  double veinAngle[5];
	
  double arteryStartPosList[4][3];
  double veinStartPosList[5][3];
	
  double arteryStartDirList[4][3];
  double veinStartDirList[5][3];
	
  // set angles eminating from backCenter, stay a distance away from edge of breast
  // edit, was 12
  double vesselEdgeSep = 20;  // distance from edge of breast of vessel entry point (mm)
	
  if(leftSide){	// left breast
    arteryAngle[0] = -4*pi/10;	// internal thoracic A
    arteryAngle[1] = 4*pi/10;	// internal thoracic B
    arteryAngle[2] = 0.0;	// thoraco-cranial
    arteryAngle[3] = -4*pi/5;	// lateral thoracic
    veinAngle[0] = pi;	// auxilluary A
    veinAngle[1] = -2*pi/3;	// auxilluary B
    veinAngle[2] = -pi/3;	// subclavian
    veinAngle[3] = pi/5;	// intercostal
    veinAngle[4] = 2*pi/3;	// internal thoracic
  } else {
    // right breast
    arteryAngle[0] = -6*pi/10;	// internal thoracic A
    arteryAngle[1] = 6*pi/10;	// internal thoracic B
    arteryAngle[2] = pi;	// thoraco-cranial
    arteryAngle[3] = -pi/5;	// lateral thoracic
    veinAngle[0] = 0.0;		// auxilluary A
    veinAngle[1] = -pi/3;		// auxilluary B
    veinAngle[2] = -2*pi/3;		// subclavian
    veinAngle[3] = 4*pi/5;		// intercostal
    veinAngle[4] = pi/3;		// internal thoracic
  }
	
  // set start positions/directions arteries
  for(int i=0; i<4; i++){
    double step = imgRes/2.0;	// step for finding boundary of breast
    double len = 0.0;
    double testPos[3];
    int voxelPos[3];
    double lcoords[3];
    testPos[0] = backCenter[0];
    bool edge = false;
    while(!edge){
      len += step;
      testPos[1] = backCenter[1] + len*cos(arteryAngle[i]);
      testPos[2] = backCenter[2] + len*sin(arteryAngle[i]);
      breast->ComputeStructuredCoordinates(testPos,voxelPos,lcoords);
      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(voxelPos));
      if(p[0] == tissue.bg || p[0] == tissue.skin){
	edge = true;
      }
    }
    arteryStartPosList[i][0] = testPos[0];
    arteryStartPosList[i][1] = backCenter[1] + (len-vesselEdgeSep)*cos(arteryAngle[i]);
    arteryStartPosList[i][2] = backCenter[2] + (len-vesselEdgeSep)*sin(arteryAngle[i]);
		
    double veclen = 0.0;
    for(int j=0; j<3; j++){
      arteryStartDirList[i][j] = nipplePos[j] - arteryStartPosList[i][j];
      veclen += arteryStartDirList[i][j]*arteryStartDirList[i][j];
    }
    veclen = sqrt(veclen);
    // normalize
    for(int j=0; j<3; j++){
      arteryStartDirList[i][j] = arteryStartDirList[i][j]/veclen;
    }
  }

  // set start positions/directions veins
  for(int i=0; i<5; i++){
    double step = imgRes/2.0;	// step for finding boundary of breast
    double len = 0.0;
    double testPos[3];
    int voxelPos[3];
    double lcoords[3];
    testPos[0] = backCenter[0];
    bool edge = false;
    while(!edge){
      len += step;
      testPos[1] = backCenter[1] + len*cos(veinAngle[i]);
      testPos[2] = backCenter[2] + len*sin(veinAngle[i]);
      breast->ComputeStructuredCoordinates(testPos,voxelPos,lcoords);
      unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(voxelPos));
      if(p[0] == tissue.bg || p[0] == tissue.skin){
	edge = true;
      }
    }
    veinStartPosList[i][0] = testPos[0];
    veinStartPosList[i][1] = backCenter[1] + (len-vesselEdgeSep)*cos(veinAngle[i]);
    veinStartPosList[i][2] = backCenter[2] + (len-vesselEdgeSep)*sin(veinAngle[i]);
		
    double veclen = 0.0;
    for(int j=0; j<3; j++){
      veinStartDirList[i][j] = nipplePos[j] - veinStartPosList[i][j];
      veclen += veinStartDirList[i][j]*veinStartDirList[i][j];
    }
    veclen = sqrt(veclen);
    // normalize
    for(int j=0; j<3; j++){
      veinStartDirList[i][j] = veinStartDirList[i][j]/veclen;
    }
  }

  // create arteries
  for(int i=0; i<4; i++){
    // create a seed for artery random number generator
    int arterySeed = static_cast<int>(round(rgen->GetRangeValue(0.0, 1.0)*2147483648));
    rgen->Next();
    
    if(i == 0){
      generate_artery(breast, vm, internalExtentVox, &tissue,
		      arteryStartPosList[i], arteryStartDirList[i], nipplePos, arterySeed, randSeed, true);
    } else {
      generate_artery(breast, vm, internalExtentVox, &tissue,
		      arteryStartPosList[i], arteryStartDirList[i], nipplePos, arterySeed, randSeed, false);
    }
  }
	
  // create veins
  for(int i=0; i<5; i++){
    // create a seed for vein random number generator
    int veinSeed = static_cast<int>(round(rgen->GetRangeValue(0.0, 1.0)*2147483648));
    rgen->Next();
		
    if(i == 0){
      generate_vein(breast, vm, internalExtentVox, &tissue,
		    veinStartPosList[i], veinStartDirList[i], nipplePos, veinSeed, randSeed, true);
    } else {
      generate_vein(breast, vm, internalExtentVox, &tissue,
		    veinStartPosList[i], veinStartDirList[i], nipplePos, veinSeed, randSeed, false);
    }
  }

  //cout << "Finished vasculature.\n";

  /*************
   * Save stuff
   ************/

  // timing
  //time(&currentTime);
  //cout << "Time interval 10: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime;


  // debug
  //time(&currentTime);
  //cout << "Time interval 20: " << difftime(currentTime,previousTime) << " seconds.\n";
  //previousTime = currentTime; 

  // save segmented breast with duct network
  vtkSmartPointer<vtkXMLImageDataWriter> writerSeg5 =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();

  writerSeg5->SetFileName(outVTIFilename);
#if VTK_MAJOR_VERSION <= 5
  writerSeg5->SetInput(breast);
#else
  writerSeg5->SetInputData(breast);
#endif
  writerSeg5->Write();

  // save metaimage format
  /*
  vtkSmartPointer<vtkMetaImageWriter> metaWrite =
    vtkSmartPointer<vtkMetaImageWriter>::New();

  metaWrite->SetCompression(1);
  metaWrite->SetRAWFileName(outrawFilename);
  metaWrite->SetFileName(outhdrFilename);

  metaWrite->SetInputData(breast);

  metaWrite->Write();
  */

  // save metaimage header
  FILE *hdrFile = fopen(outhdrFilename, "w");

  if(hdrFile == NULL){
    cerr << "Unable to open mhd file for writing\n";
  } else {
    fprintf(hdrFile, "ObjectType = Image\n");
    fprintf(hdrFile, "NDims = 3\n");
    fprintf(hdrFile, "BinaryData = True\n");
    fprintf(hdrFile, "BinaryDataByteOrderMSB = False\n");
    fprintf(hdrFile, "CompressedData = False\n");
    fprintf(hdrFile, "TransformMatrix = 1 0 0 0 1 0 0 0 1\n");
    fprintf(hdrFile, "Offset = %6.4f %6.4f %6.4f\n", origin[0], origin[1], origin[2]);
    fprintf(hdrFile, "CenterOfRotation = 0 0 0\n");
    fprintf(hdrFile, "ElementSpacing = %6.4f %6.4f %6.4f\n", spacing[0], spacing[1], spacing[2]);
    fprintf(hdrFile, "DimSize = %d %d %d\n", dim[0], dim[1], dim[2]);
    fprintf(hdrFile, "AnatomicalOrientation = ???\n");
    fprintf(hdrFile, "ElementType = MET_UCHAR\n");
    fprintf(hdrFile, "ObjectType = Image\n");
    fprintf(hdrFile, "ElementDataFile = p_%d.raw\n", randSeed);
    fclose(hdrFile);
  }

  // save gzipped raw
  gzFile gzf = gzopen(outgzFilename, "wb");

  if(gzf == NULL){
    cerr << "Unable to open gzip file for writing\n";
  } else {
    gzbuffer(gzf, dim[1]*dim[2]);
    unsigned char *p = static_cast<unsigned char*>(breast->GetScalarPointer());
    for(int i=0; i<dim[0]; i++){
      gzwrite(gzf, static_cast<const void*>(p), dim[1]*dim[2]);
      p += dim[1]*dim[2];
    }

    gzclose(gzf);
  }




  // save some projections in raw format

  /*

  unsigned int* skinProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* fatProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* ligProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* ductProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* glandProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* TDLUProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* arteryProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* veinProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));
  unsigned int* nippleProj = (unsigned int*)calloc(dim[0]*dim[2],sizeof(unsigned int));

  for(int i=0; i<dim[0]; i++){
    for(int k=0; k<dim[2]; k++){
      for(int j=0; j<dim[1]; j++){
	unsigned char* p = static_cast<unsigned char*>(breast->GetScalarPointer(i,j,k));
	if(p[0] == tissue.cooper){
	  ligProj[i*dim[2]+k] = ligProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.fat){
	  fatProj[i*dim[2]+k] = fatProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.skin){
	  skinProj[i*dim[2]+k] = skinProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.duct){
	  ductProj[i*dim[2]+k] = ductProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.TDLU){
	  TDLUProj[i*dim[2]+k] = TDLUProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.gland){
	  glandProj[i*dim[2]+k] = glandProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.artery){
	  arteryProj[i*dim[2]+k] = arteryProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.vein){
	  veinProj[i*dim[2]+k] = veinProj[i*dim[2]+k] + 1;
	} else if(p[0] == tissue.nipple){
	  nippleProj[i*dim[2]+k] = nippleProj[i*dim[2]+k] + 1;
	}
      }
    }
  }

  cout << "Output projection dimensions = " << dim[0] << " by " << dim[2] << "\n";
    
  // save in raw format

  FILE *projFile;

  char skinProjFilename[128];
  sprintf(skinProjFilename,"%s/p_%d_%d_%d_skin.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char fatProjFilename[128];
  sprintf(fatProjFilename,"%s/p_%d_%d_%d_fat.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char ligProjFilename[128];
  sprintf(ligProjFilename,"%s/p_%d_%d_%d_lig.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char ductProjFilename[128];
  sprintf(ductProjFilename,"%s/p_%d_%d_%d_duct.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char glandProjFilename[128];
  sprintf(glandProjFilename,"%s/p_%d_%d_%d_gland.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char TDLUProjFilename[128];
  sprintf(TDLUProjFilename,"%s/p_%d_%d_%d_TDLU.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char arteryProjFilename[128];
  sprintf(arteryProjFilename,"%s/p_%d_%d_%d_artery.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char veinProjFilename[128];
  sprintf(veinProjFilename,"%s/p_%d_%d_%d_vein.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);
  char nippleProjFilename[128];
  sprintf(nippleProjFilename,"%s/p_%d_%d_%d_nipple.dat", outputDir.c_str(),randSeed,dim[0],dim[2]);

  projFile = fopen(skinProjFilename,"wb");
  fwrite(skinProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(fatProjFilename,"wb");
  fwrite(fatProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(ligProjFilename,"wb");
  fwrite(ligProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(ductProjFilename,"wb");
  fwrite(ductProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(glandProjFilename,"wb");
  fwrite(glandProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(TDLUProjFilename,"wb");
  fwrite(TDLUProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(arteryProjFilename,"wb");
  fwrite(arteryProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(veinProjFilename,"wb");
  fwrite(veinProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);
  projFile = fopen(nippleProjFilename,"wb");
  fwrite(nippleProj, dim[0]*dim[2], sizeof(unsigned int), projFile);
  fclose(projFile);

  free(skinProj);
  free(fatProj);
  free(ligProj);
  free(ductProj);
  free(glandProj);
  free(TDLUProj);
  free(arteryProj);
  free(veinProj);
  free(nippleProj);
  */
  // save actual vessel distance maps

/*
	
  vtkSmartPointer<vtkXMLImageDataReader> fillReader =
  vtkSmartPointer<vtkXMLImageDataReader>::New();
  fillReader->SetFileName("veinFill.vti");
  fillReader->Update();
  vtkSmartPointer<vtkImageData> holder = 
  vtkSmartPointer<vtkImageData>::New();
  holder = fillReader->GetOutput();
  
  double* fd = static_cast<double*>(holder->GetScalarPointer());
  for(int i=0; i<holder->GetNumberOfPoints(); i++){
  fd[i] = sqrt(fd[i]);
  }
	
  vtkSmartPointer<vtkXMLImageDataWriter> fillWriter =
  vtkSmartPointer<vtkXMLImageDataWriter>::New();
  
  fillWriter->SetFileName("veinDist.vti");
  #if VTK_MAJOR_VERSION <= 5
  fillWriter->SetInput(holder);
  #else
  fillWriter->SetInputData(holder);
  #endif
  fillWriter->Write();
	
  fillReader->SetFileName("arteryFill.vti");
  fillReader->Update();
  holder = fillReader->GetOutput();
	
  fd = static_cast<double*>(holder->GetScalarPointer());

  for(int i=0; i<holder->GetNumberOfPoints(); i++){
  fd[i] = sqrt(fd[i]);
  }
	
  fillWriter->SetFileName("arteryDist.vti");
  #if VTK_MAJOR_VERSION <= 5
  fillWriter->SetInput(holder);
  #else
  fillWriter->SetInputData(holder);
  #endif
  fillWriter->Write();

  // debug
  vtkSmartPointer<vtkTIFFWriter> tiffWriter =
  vtkSmartPointer<vtkTIFFWriter>::New();
  tiffWriter->SetCompressionToPackBits();
  tiffWriter->SetFileName("final.tiff");
  #if VTK_MAJOR_VERSION <= 5
  tiffWriter->SetInput(breast);
  #else
  tiffWriter->SetInputData(breast);
  #endif
  tiffWriter->Write();
*/
	

//time(&endTime);
//	cout << "Phantom generation time: " << difftime(endTime,startTime) << " seconds\n";

  return EXIT_SUCCESS;
}

inline void statusBar(unsigned int current, unsigned int total, unsigned int width, unsigned int numUpdate){
	
  double prevNum = floor((double)(current-1)/(double)(total)*(double)(numUpdate));
  double curNum = floor((double)(current)/(double)(total)*(double)(numUpdate));

  if(curNum != prevNum){
    // draw new status bar
    cout << setw(3) << (int)((double)(current)/(double)(total)*100) << "% [";
    for(int x=0; x<=(int)curNum; x++){
      cout << "=";
    }
    for(int x=(int)curNum+1; x<=width; x++){
      cout << " ";
    }
    cout << "]\r" << std::flush;
  }
}


