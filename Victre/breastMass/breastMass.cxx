/*! \file breastMass.cxx
 *  \brief breastMass main source file
 *  \author Christian G. Graff
 *  \version 1.0
 *  \date 2018
 *
 * based on "A computational model to generate simulated three-dimensional breast masses"
 * Sisternes et al., Med. Phys 42(2) 2015.
 *  
 *  \copyright To the extent possible under law, the author(s) have
 *  dedicated all copyright and related and neighboring rights to this
 *  software to the public domain worldwide. This software is
 *  distributed without any warranty.  You should have received a copy
 *  of the CC0 Public Domain Dedication along with this software.
 *  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 * 
 */

#include "breastMass.hxx"

namespace po = boost::program_options;

struct perturb {
  double r;
  double theta;
  double phi;
  double length;
  double rad;
  bool spike;
};

int main(int argc, char* argv[]){

  double pi = vtkMath::Pi();
	
  // configuration variables
  po::options_description baseOpt("base options");
  baseOpt.add_options()
    ("base.imgRes",po::value<double>()->default_value(0.25),"voxel size (mm)")
    ("base.seed",po::value<unsigned int>(),"random number generator seed")
    ("base.complexity",po::value<double>()->default_value(1.0), "complexity scaling, 1=normal")
    ;
	
  po::options_description massOpt("mass options");
  massOpt.add_options()
    ("mass.lMax",po::value<int>()->default_value(3),"maximum spherical harmonic order")
    ("mass.alpha", po::value<double>()->default_value(10.0),"mean mass radius (mm)")
    ("mass.meanSigma2", po::value<double>()->default_value(0.31),"mean mass surface irregularity variance (mm^2)")
    ("mass.stdSigma2", po::value<double>()->default_value(0.04),"mass surface irregularity standard deviation (mm^2)")
    ("mass.powerLaw",po::value<double>()->default_value(4.0),"covariance power law index")
    ("mass.meanLF",po::value<double>()->default_value(611.2),"mean number of low freq. modifications")
    ("mass.stdLF",po::value<double>()->default_value(70.6),"std. deviation of number of low freq. modifications")
    ("mass.meanShape",po::value<double>()->default_value(0.36),"mean LF shape distribution 0=spike, 1=bump")
    ("mass.stdShape",po::value<double>()->default_value(0.48),"std. deviation of LF shape 0=spike, 1=bump")
    ("mass.meanLFRad",po::value<double>()->default_value(0.229),"mean relative LF radius")
    ("mass.stdLFRad",po::value<double>()->default_value(0.073),"std. deviation of LF radius")
    ("mass.meanLFLen",po::value<double>()->default_value(0.113),"mean relative LF length")
    ("mass.stdLFLen",po::value<double>()->default_value(0.021),"std. deviation of LF length")
    ("mass.meanFuzzAlpha",po::value<double>()->default_value(0.015),"mean fuzzy alpha")
    ("mass.stdFuzzAlpha",po::value<double>()->default_value(0.0),"std. deviation of fuzzy alpha")
    ;

  po::options_description spiculeOpt("spicule options");
  spiculeOpt.add_options()
    ("spicule.meanInitial",po::value<double>()->default_value(0),"mean number of initial segments")
    ("spicule.stdInitial",po::value<double>()->default_value(0),"std. deviation number of initial segments")
    ("spicule.meanNeigh",po::value<double>()->default_value(8.98),"mean max number of neighbor segments")
    ("spicule.stdNeigh",po::value<double>()->default_value(1.89),"std. deviation max number of neighbor segments")
    ("spicule.meanInitRad",po::value<double>()->default_value(0.024),"mean initial relative radius")
    ("spicule.stdInitRad",po::value<double>()->default_value(0.0053),"std. deviation initial relative radius")
    ("spicule.meanRadDec",po::value<double>()->default_value(0.89),"mean radius decrease")
    ("spicule.stdRadDec",po::value<double>()->default_value(0.31),"std. deviation radius decrease")
    ("spicule.meanInitLen",po::value<double>()->default_value(0.173),"mean initial relative radius")
    ("spicule.stdInitLen",po::value<double>()->default_value(0.018),"std. deviation initial relative radius")
    ("spicule.meanLenDec",po::value<double>()->default_value(0.91),"mean length decrease")
    ("spicule.stdLenDec",po::value<double>()->default_value(0.30),"std. deviation length decrease")
    ("spicule.meanContProb",po::value<double>()->default_value(0.717),"mean continue prob.")
    ("spicule.stdContProb",po::value<double>()->default_value(0.057),"std. deviation continue prob.")
    // Assuming equi-probably symmetric and asymmetric branching
    //("spicule.meanSymBifProb",po::value<double>()->default_value(0.142),"mean symmentric bifurcation prob.")
    //("spicule.stdSymBifProb",po::value<double>()->default_value(0.028),"std. deviation symm bif. prob.")
    //("spicule.meanAsymBifProb",po::value<double>()->default_value(0.142),"mean asymmentric bifurcation prob.")
    //("spicule.stdAsymBifProb",po::value<double>()->default_value(0.028),"std. deviation asymm bif. prob.")
    ("spicule.meanBranchAng",po::value<double>()->default_value(6.55),"mean asymmentric bifurcation prob.")
    ("spicule.stdBranchAng",po::value<double>()->default_value(0.62),"std. deviation asymm bif. prob.")
    ;
		
  // config file options
  po::options_description configFileOpt("config file options");
  configFileOpt.add(baseOpt);
  configFileOpt.add(massOpt);
  configFileOpt.add(spiculeOpt);

  // all of the options
  po::options_description all("All options");
  all.add_options()
    ("config,c", po::value<std::string>()->required(), "name of configuration file")
    ("help,h", "print help message")
    ;
  all.add(configFileOpt);

  po::variables_map vm;

  // get configuration filename from command line
  po::store(parse_command_line(argc,argv,all), vm);

  if (vm.count("help")) {
    cout << "\n";
    cout << "breastMass - breast mass generation\n\n";
    // cout << opts << "\n";
    return 1;
  }
  
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
    size_t result;
    result = fread((char*)(&randSeed), sizeof(int),1,randFile);
    fclose(randFile);
  }
	
  // start a random number generator
  vtkSmartPointer<vtkMinimalStandardRandomSequence> rgen =
    vtkSmartPointer<vtkMinimalStandardRandomSequence>::New();

  rgen->SetSeed(randSeed);
	
  // start a standard normal random number generator
  vtkSmartPointer<vtkBoxMuellerRandomSequence> rNormGen = 
    vtkSmartPointer<vtkBoxMuellerRandomSequence>::New();
		
  rNormGen->SetUniformSequence(rgen);

  double imgRes = vm["base.imgRes"].as<double>();
  double complexity = vm["base.complexity"].as<double>();

  // setup GRS parameters
	
  double alpha = vm["mass.alpha"].as<double>();
  int lMax = vm["mass.lMax"].as<int>();
  double nu = vm["mass.powerLaw"].as<double>();
  double sigma2 = rNormGen->GetScaledValue(vm["mass.meanSigma2"].as<double>()*complexity,
					   vm["mass.stdSigma2"].as<double>());
	
  if(sigma2 < 0.0){
    sigma2 = 0.0;
  }
	
  double pLawSum = 0.0;
	
  // note: error in paper - index start
  for(int i=1; i<=lMax; i++){
    pLawSum += 1.0/pow(i, nu);
  }
	
  double Ctilde = log(1+sigma2)*(1.0/pLawSum);
	
  // random spherical harmonic coefficients
  std::complex<double>** slm = new std::complex<double>*[lMax+1];
  for(int i=0; i<=lMax; i++){
    slm[i] = new std::complex<double>[lMax+1];
  }
		
  for(int l=0; l<=lMax; l++){
    double Cl;
    if(l == 0 || l == 1){
      Cl = 0.0;
    } else {
      Cl = Ctilde/pow(l, nu);
    }	
    for(int m=0; m<=l; m++){
      if(m == 0){
	double tempValr = rNormGen->GetValue()*sqrt(2.0);
	rNormGen->Next();
	slm[l][m] = std::complex<double>(tempValr, 0.0);
      } else {
	double tempValr = rNormGen->GetValue();
	rNormGen->Next();
	double tempVali = rNormGen->GetValue();
	rNormGen->Next();
	slm[l][m] = std::complex<double>(tempValr, tempVali);
      }
      
      slm[l][m] = slm[l][m]*sqrt(2*pi*Cl/(2*l+1));
    }
  }
	
  // determine spikes and bumps
  int numLF = static_cast<int>(round(rNormGen->GetScaledValue(vm["mass.meanLF"].as<double>(), 
							      vm["mass.stdLF"].as<double>())));
  rNormGen->Next();
  
  if(numLF < 0){
    numLF = 0;
  }
	
  double typeFrac = rNormGen->GetScaledValue(vm["mass.meanShape"].as<double>(),
					     vm["mass.stdShape"].as<double>());
  rNormGen->Next();
  
  if(typeFrac < 0.0){
    typeFrac = 0.0;
  }
  if(typeFrac > 1.0){
    typeFrac = 1.0;
  }
	
  int numBump = static_cast<int>(round(numLF*typeFrac));
  int numSpike = numLF - numBump;
	
  // pick locations, sizes
  perturb *LF = (perturb*)malloc(numLF*sizeof(perturb));
	
  for(int i=0; i<numLF; i++){
    LF[i].length = rNormGen->GetScaledValue(vm["mass.meanLFLen"].as<double>()*complexity,
					    vm["mass.stdLFLen"].as<double>())*alpha;
    rNormGen->Next();
    if(LF[i].length < 0.0){
      LF[i].length = 0.0;
    }
    LF[i].rad = rNormGen->GetScaledValue(vm["mass.meanLFRad"].as<double>(),
					 vm["mass.stdLFRad"].as<double>())*alpha;
    rNormGen->Next();
    
    LF[i].spike = (i < numSpike) ? true : false; 
		
    // random point on sphere
    LF[i].phi = rgen->GetRangeValue(-pi, pi);
    rgen->Next();
    LF[i].theta = acos(2*rgen->GetValue()-1);
    rgen->Next();
		
    // GRS radius
    std::complex<double> s(0.0, 0.0);
    for(int l=0; l<=lMax; l++){
      for(int m=0; m<=l; m++){
	s = s + slm[l][m] * boost::math::spherical_harmonic(l, m, LF[i].theta, LF[i].phi);
				
	if(m > 0){
	  // do negative index
	  s = s + pow(-1.0, m)*std::conj(slm[l][m]) * boost::math::spherical_harmonic(l, -m, LF[i].theta, LF[i].phi);
	}
      }
    }
				
    std::complex<double> sexp;
    sexp = std::exp(s);
				
    LF[i].r = alpha/sqrt(1.0+sigma2)*sexp.real();
  }

  // create radial extent map
	
  int ustep = (int)ceil(alpha*2*pi/imgRes);
  int vstep = (int)ceil(alpha*pi/imgRes);
	
  // track max value
  double rmax = 0.0;
	
  vtkSmartPointer<vtkImageData> rmap = 
    vtkSmartPointer<vtkImageData>::New();
	
  double rmapSpacing[3] = {1.0/double(ustep-1), 1.0/double(vstep-1), 1};
	
  rmap->SetSpacing(rmapSpacing);
  rmap->SetExtent(0, ustep-1, 0, vstep-1, 0, 0);
  rmap->SetOrigin(0.0,0.0,0.0);
#if VTK_MAJOR_VERSION <= 5
  rmap->SetNumberOfScalarComponents(1);
  rmap->SetScalarTypeToDouble();
  rmap->AllocateScalars();
#else
  rmap->AllocateScalars(VTK_DOUBLE,1);
#endif
	
  // fill in rmap
  for(int i=0; i<ustep; i++){
    double uval = 0+i*rmapSpacing[0];
    for(int j=0; j<vstep; j++){
      double vval = 0+j*rmapSpacing[1];
      
      double theta = acos(2*vval-1.0);
      double phi = -pi + 2*pi*uval;
			
      // GRS model
      std::complex<double> s(0.0, 0.0);
      for(int l=0; l<=lMax; l++){
	for(int m=0; m<=l; m++){
	  s = s + slm[l][m] * boost::math::spherical_harmonic(l, m, theta, phi);
					
	  if(m > 0){
	    // do negative index
	    s = s + pow(-1.0, m)*std::conj(slm[l][m]) * boost::math::spherical_harmonic(l, -m, theta, phi);
	  }
	}
      }
			
      std::complex<double> sexp;
      sexp = std::exp(s);
			
      // GRS model radius
      double rGRS = alpha/sqrt(1.0+sigma2)*sexp.real();
			
      // check for LF perturbation
      double rLF = rGRS;
			
#pragma omp parallel for
      for(int a=0; a<numLF; a++){
	double dr;
				
	// check distance
	double dist = pow(rGRS*sin(theta)*sin(phi) - LF[a].r*sin(LF[a].theta)*sin(LF[a].phi),2.0);
	dist += pow(rGRS*sin(theta)*cos(phi) - LF[a].r*sin(LF[a].theta)*cos(LF[a].phi),2.0);
	dist += pow(rGRS*cos(theta) - LF[a].r*cos(LF[a].theta),2.0);
				
	dist = sqrt(dist);
	if(dist <= LF[a].rad){
	  // perturbation
	  if(LF[a].spike){
	    dr = LF[a].length*pow((LF[a].rad - dist)/LF[a].rad, 4.0);
	  } else {
	    // bump
	    dr = LF[a].length*(exp(-dist/2.0/pow(LF[a].rad/2.0,2)) - exp(-4.5));
	  }
					
#pragma omp atomic
	  rLF += dr;
	}
      }
			
      double rFuzz = rLF;
      
      double fuzzAlpha = rNormGen->GetScaledValue(vm["mass.meanFuzzAlpha"].as<double>(),
						  vm["mass.stdFuzzAlpha"].as<double>());
      rNormGen->Next();
			
      rFuzz = rFuzz*(1.0 + fuzzAlpha*rNormGen->GetValue());
      rNormGen->Next();
			
      // save in rmap
      double* p = static_cast<double*>(rmap->GetScalarPointer(i,j,0));
			
      *p = rFuzz;
			
      if(rFuzz > rmax){
	rmax = rFuzz;
      }
    }
  }

  // create mass
	
  // define mass data structure
  vtkSmartPointer<vtkImageData> mass = 
    vtkSmartPointer<vtkImageData>::New();
	
  double massSpacing[3] = {imgRes, imgRes, imgRes};
	
  mass->SetSpacing(massSpacing);
	
  // set mass extent based on rmax
  int massDim = static_cast<int>(ceil(2*2*rmax/imgRes));
	
	
  mass->SetExtent(-massDim/2, -massDim/2+massDim-1, -massDim/2, -massDim/2+massDim-1, -massDim/2, -massDim/2+massDim-1);
  mass->SetOrigin(0.0,0.0,0.0);
	
#if VTK_MAJOR_VERSION <= 5
  mass->SetNumberOfScalarComponents(1);
  mass->SetScalarTypeToUnsignedChar();
  mass->AllocateScalars();
#else
  mass->AllocateScalars(VTK_UNSIGNED_CHAR,1);
#endif
	
  int massExtent[6];
  mass->GetExtent(massExtent);
	
  vtkSmartPointer<vtkImageInterpolator> inter = 
    vtkSmartPointer<vtkImageInterpolator>::New();
		
  inter->Initialize(rmap);
  inter->SetInterpolationModeToLinear();
	
#pragma omp parallel for
  for(int i=massExtent[0]; i<=massExtent[1]; i++){
    for(int j=massExtent[2]; j<=massExtent[3]; j++){
      for(int k=massExtent[4]; k<=massExtent[5]; k++){
	// spherical coordinates
	double r = sqrt(i*i+j*j+k*k)*imgRes;
	double theta = acos(k*imgRes/r);
	double phi = atan2(j*imgRes, i*imgRes);
				
	double uval = 0.5+phi/2.0/pi;
	double vval = (cos(theta)+1.0)/2.0;
				
	// interpolate rmap
	double rval;
	double coord[3] = {uval, vval, 0.0};
	inter->Interpolate(coord, &rval);
				
	unsigned char* p = static_cast<unsigned char*>(mass->GetScalarPointer(i,j,k));
	// test if inside mass
	if(r <= rval){
	  // in mass
	  *p = 1;
	} else {
	  *p = 0;
	}
      }
    }
  }
		
  // spiculations
	
  int numInitial = static_cast<int>(round(rNormGen->GetScaledValue(vm["spicule.meanInitial"].as<double>(),
								   vm["spicule.stdInitial"].as<double>())));
	
	
  // spicule continuation probability
  double contProb = rNormGen->GetScaledValue(vm["spicule.meanContProb"].as<double>(),
					     vm["spicule.stdContProb"].as<double>());
	
  if(contProb > 1.0){
    contProb = 1.0;
  }
  if(contProb < 0.0){
    contProb = 0.0;
  }
	
  // bifurcation angle variance
  double gamma = rNormGen->GetScaledValue(vm["spicule.meanBranchAng"].as<double>(),
				vm["spicule.stdBranchAng"].as<double>())*pi/180.0;
  
  if(gamma > pi){
    gamma = pi;
  }
  if(gamma < 0.0){
    gamma = 0.0;
  }
	
  int initialSeg = 0;

  while(initialSeg<numInitial){
    // start new neighborhood
		
    // location
    double phiIni = rgen->GetRangeValue(-pi, pi);
    rgen->Next();
    double thetaIni = acos(2*rgen->GetValue()-1);
    rgen->Next();
		
    // number in neighborhood
    int numSeg = static_cast<int>(round(rNormGen->GetScaledValue(vm["spicule.meanNeigh"].as<double>(),
								 vm["spicule.stdNeigh"].as<double>())));
				
    if(numSeg<0){
      numSeg = 0;
    }
		
    int nDone = 0;
		
    while(nDone<numSeg && initialSeg<numInitial){
      // branch angle
			
      double phi = phiIni + (floor(log2(nDone+1))+rgen->GetValue())*pi/50.0;
      rgen->Next();
      double theta = thetaIni + (floor(log2(nDone+1))+rgen->GetValue())*pi/50.0;
      rgen->Next();
			
      double rad = rNormGen->GetScaledValue(vm["spicule.meanInitRad"].as<double>(),
					    vm["spicule.stdInitRad"].as<double>())*alpha;
      rNormGen->Next();
				
      if(rad<0.0){
	rad=0.0;
      }
			
      double len = rNormGen->GetScaledValue(vm["spicule.meanInitLen"].as<double>(),
					    vm["spicule.stdInitLen"].as<double>())*alpha;
      rNormGen->Next();
				
      if(len<0.0){
	len=0.0;
      }
		
      double thetaDir = theta + rNormGen->GetScaledValue(0.0, 15.0*pi/180.0);
      rNormGen->Next();
      double phiDir = phi + rNormGen->GetScaledValue(0.0, 15.0*pi/180.0);
      rNormGen->Next();
		
      // radial position
			
      double uval = 0.5+phi/2.0/pi;
      double vval = (cos(theta)+1.0)/2.0;
				
      // interpolate rmap
      double r;
      double coord[3] = {uval, vval, 0.0};
      inter->Interpolate(coord, &r);
			
      createBranch(r, theta, phi, thetaDir, phiDir, rad, len, contProb, gamma, mass, rgen, rNormGen, vm);
			
      nDone += 1;
      initialSeg += 1;
    }
  }

	
  // outputs
	
  // filename stub
  char fileStub[128];
	
  sprintf(fileStub, "mass_%i", randSeed);
	
	
  // copy configuration
  char cfgFilename[128];
  sprintf(cfgFilename, "%s.cfg", fileStub);
  FILE *cfgCopy = fopen(cfgFilename,"wt");
  FILE *cfgRead = fopen(configFile.c_str(),"rt");
	
  int readChar = getc(cfgRead);
  while(readChar != EOF){
    putc(readChar, cfgCopy);
    readChar = getc(cfgRead);
  }
  fclose(cfgCopy);
  fclose(cfgRead);

  // save volume
  vtkSmartPointer<vtkXMLImageDataWriter> writer =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();

  char volFilename[128];
  sprintf(volFilename, "%s.vti", fileStub);
  writer->SetFileName(volFilename);
#if VTK_MAJOR_VERSION <= 5
  writer->SetInput(mass);
#else
  writer->SetInputData(mass);
#endif
  writer->Write();
		
  // save volume raw format
  char rawFilename[128];
  sprintf(rawFilename, "%s_%i.raw", fileStub, massDim);
	
  FILE *rawOut;
	
  rawOut = fopen(rawFilename,"wb");
	
  unsigned char* out = static_cast<unsigned char*>(mass->GetScalarPointer());
  fwrite(out, massDim*massDim*massDim, 1, rawOut);
  fclose(rawOut);
	
  // free slm
  for(int i=0; i<=lMax; i++){
    delete[] slm[i];
  }
  delete[] slm;
	
  return EXIT_SUCCESS;
}
					
void createBranch(double r, double theta, double phi, double thetaDir, double phiDir, 
		  double rad, double len, double contProb, double gam, vtkImageData* mass, 
		  vtkMinimalStandardRandomSequence* rgen, vtkBoxMuellerRandomSequence* rNormGen, po::variables_map vm){

  double pi = vtkMath::Pi();
  double imgRes = vm["base.imgRes"].as<double>();

  // end radius
  double endRad = rNormGen->GetScaledValue(vm["spicule.meanRadDec"].as<double>(),
					   vm["spicule.stdRadDec"].as<double>())*rad;
	
  if(endRad > rad){
    endRad = rad;
  }
	
  if(endRad < 0.0){
    endRad = 0.0;
  }
	
  // segment starting sphere
  int startPix[3];
  double startPos[3] = {r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta)};
  double pcoords[3];
  int inROI;
  int searchBox[6];
	
  inROI = mass->ComputeStructuredCoordinates(startPos, startPix, pcoords);
	
  if(inROI){
		
    //  fill starting sphere
    int searchRad = static_cast<int>(ceil(rad/imgRes));
    int endRadPix = static_cast<int>(ceil(endRad/imgRes));
		
    int massExtent[6];
    mass->GetExtent(massExtent);
		
    searchBox[0] = (startPix[0]-searchRad < massExtent[0]) ? massExtent[0] :  startPix[0]-searchRad;
    searchBox[1] = (startPix[0]+searchRad > massExtent[1]) ? massExtent[1] :  startPix[0]+searchRad;
    searchBox[2] = (startPix[1]-searchRad < massExtent[2]) ? massExtent[2] :  startPix[1]-searchRad;
    searchBox[3] = (startPix[1]+searchRad > massExtent[3]) ? massExtent[3] :  startPix[1]+searchRad;
    searchBox[4] = (startPix[2]-searchRad < massExtent[4]) ? massExtent[4] :  startPix[2]-searchRad;
    searchBox[5] = (startPix[2]+searchRad > massExtent[5]) ? massExtent[5] :  startPix[2]+searchRad;
		
    for(int a=searchBox[0]; a<=searchBox[1]; a++){
      for(int b=searchBox[2]; b<=searchBox[3]; b++){
	for(int c=searchBox[4]; c<=searchBox[5]; c++){
	  double dist = imgRes*sqrt((a-startPix[0])*(a-startPix[0]) + (b-startPix[1])*(b-startPix[1]) + (c-startPix[2])*(c-startPix[2]));
	  if(dist <= rad){
	    unsigned char* p = static_cast<unsigned char*>(mass->GetScalarPointer(a,b,c));
	    *p = 1;
	  }
	}
      }
    }
		
    // find end point
    double endPos[3];
    endPos[0] = startPos[0] + len*sin(thetaDir)*cos(phiDir);
    endPos[1] = startPos[1] + len*sin(thetaDir)*sin(phiDir);
    endPos[2] = startPos[2] + len*cos(thetaDir);
    
    int endPix[3];
    int endOK;
    endOK = mass->ComputeStructuredCoordinates(endPos, endPix, pcoords);
		
    // unit vector
    double uvec[3] = {sin(thetaDir)*cos(phiDir), sin(thetaDir)*sin(phiDir), cos(thetaDir)};
		
    // determine search box for segment
    searchBox[0] = (searchBox[0] < endPix[0]-endRadPix) ? searchBox[0] : endPix[0]-endRadPix;
    searchBox[1] = (searchBox[1] > endPix[0]+endRadPix) ? searchBox[1] : endPix[0]+endRadPix;
    searchBox[2] = (searchBox[2] < endPix[1]-endRadPix) ? searchBox[2] : endPix[1]-endRadPix;
    searchBox[3] = (searchBox[3] > endPix[1]+endRadPix) ? searchBox[3] : endPix[1]+endRadPix;
    searchBox[4] = (searchBox[4] < endPix[2]-endRadPix) ? searchBox[4] : endPix[2]-endRadPix;
    searchBox[5] = (searchBox[5] > endPix[2]+endRadPix) ? searchBox[5] : endPix[2]+endRadPix;
		
    searchBox[0] = (searchBox[0] < massExtent[0]) ? massExtent[0] :  searchBox[0];
    searchBox[1] = (searchBox[1] > massExtent[1]) ? massExtent[1] :  searchBox[1];
    searchBox[2] = (searchBox[2] < massExtent[2]) ? massExtent[2] :  searchBox[2];
    searchBox[3] = (searchBox[3] > massExtent[3]) ? massExtent[3] :  searchBox[3];
    searchBox[4] = (searchBox[4] < massExtent[4]) ? massExtent[4] :  searchBox[4];
    searchBox[5] = (searchBox[5] > massExtent[5]) ? massExtent[5] :  searchBox[5];
		
    // fill segment
#pragma omp parallel for
    for(int a=searchBox[0]; a<=searchBox[1]; a++){
      for(int b=searchBox[2]; b<=searchBox[3]; b++){
	for(int c=searchBox[4]; c<=searchBox[5]; c++){
	  double curPos[3] = {a*imgRes, b*imgRes, c*imgRes};
	  // project to center line
	  double ldist = uvec[0]*(curPos[0]-startPos[0])+uvec[1]*(curPos[1]-startPos[1])+uvec[2]*(curPos[2]-startPos[2]);
									
	  if(ldist <= len && ldist >= 0.0){
	    // current radius
	    double thisRad = rad+(endRad-rad)*ldist/len;
						
	    // radial distance - pythagoras
	    double curRad = sqrt((curPos[0]-startPos[0])*(curPos[0]-startPos[0])+(curPos[1]-startPos[1])*(curPos[1]-startPos[1])
				 +(curPos[2]-startPos[2])*(curPos[2]-startPos[2])-ldist*ldist);
	    
	    if(curRad <= thisRad){
	      unsigned char* p = static_cast<unsigned char*>(mass->GetScalarPointer(a,b,c));
	      *p = 1;
	    }
	  }
	}
      }
    }
		
    //  fill ending sphere
		
    searchBox[0] = (endPix[0]-endRadPix < massExtent[0]) ? massExtent[0] :  endPix[0]-endRadPix;
    searchBox[1] = (endPix[0]+endRadPix > massExtent[1]) ? massExtent[1] :  endPix[0]+endRadPix;
    searchBox[2] = (endPix[1]-endRadPix < massExtent[2]) ? massExtent[2] :  endPix[1]-endRadPix;
    searchBox[3] = (endPix[1]+endRadPix > massExtent[3]) ? massExtent[3] :  endPix[1]+endRadPix;
    searchBox[4] = (endPix[2]-endRadPix < massExtent[4]) ? massExtent[4] :  endPix[2]-endRadPix;
    searchBox[5] = (endPix[2]+endRadPix > massExtent[5]) ? massExtent[5] :  endPix[2]+endRadPix;
		
    for(int a=searchBox[0]; a<=searchBox[1]; a++){
      for(int b=searchBox[2]; b<=searchBox[3]; b++){
	for(int c=searchBox[4]; c<=searchBox[5]; c++){
	  double dist = imgRes*sqrt((a-endPix[0])*(a-endPix[0]) + (b-endPix[1])*(b-endPix[1]) + (c-endPix[2])*(c-endPix[2]));
	  if(dist <= endRad){
	    unsigned char* p = static_cast<unsigned char*>(mass->GetScalarPointer(a,b,c));
	    *p = 1;
	  }
	}
      }
    }
				
    // create children
    if(endRad>imgRes/2.0 && endOK >= 0){
      double newR = sqrt(endPos[0]*endPos[0]+endPos[1]*endPos[1]+endPos[2]*endPos[2]);
      double newTheta = acos(endPos[2]/newR);
      double newPhi = atan2(endPos[1], endPos[0]);
      double newLen = len*rNormGen->GetScaledValue(vm["spicule.meanLenDec"].as<double>(),
						   vm["spicule.stdLenDec"].as<double>());
      rNormGen->Next();
      
      if(newLen > len){
	newLen = len;
      }
			
      if(newLen < 0.0){
	newLen = 0.0;
      }
			
      double rval = rgen->GetValue();
      rgen->Next();
      
      if(rval < contProb){
	// one child
				
	// branching angle
	double newThetaDir =  thetaDir + gam*rNormGen->GetValue();
	rNormGen->Next();
	double newPhiDir =  phiDir + gam*rNormGen->GetValue();
	rNormGen->Next();
				
	// create
	createBranch(newR, newTheta, newPhi, newThetaDir, newPhiDir, endRad, newLen, contProb, gam, mass, rgen, rNormGen, vm);
				
      } else {
	// two children
	double rval2 = rgen->GetValue();
	rgen->Next();
				
	if(rval2 < 0.5){
	  // symmetric bifurcation
	  double dr = rgen->GetRangeValue(0.5, 0.8);
	  rgen->Next();
	  double newThetaDir = thetaDir + gam*(rNormGen->GetValue()+1.0);
	  rNormGen->Next();
	  double newPhiDir = phiDir + gam*rNormGen->GetValue();
	  rNormGen->Next();
	  createBranch(newR, newTheta, newPhi, newThetaDir, newPhiDir, endRad*pow(dr,1.0/2.6), newLen*pow(dr,1.0/2.6), contProb, gam, mass, rgen, rNormGen, vm);
	  newThetaDir = thetaDir + gam*(rNormGen->GetValue()-1.0);
	  rNormGen->Next();
	  newPhiDir = phiDir + gam*rNormGen->GetValue();
	  rNormGen->Next();
	  createBranch(newR, newTheta, newPhi, newThetaDir, newPhiDir, endRad*pow(1-dr,1.0/2.6), newLen*pow(1-dr,1.0/2.6), contProb, gam, mass, rgen, rNormGen, vm);
					
	} else {
	  // asymmetric bifurcation
	  double dr = rgen->GetRangeValue(0.8, 1.0);
	  rgen->Next();
	  double newThetaDir = thetaDir + gam*rNormGen->GetValue();
	  rNormGen->Next();
	  double newPhiDir = phiDir + gam*rNormGen->GetValue();
	  rNormGen->Next();
	  createBranch(newR, newTheta, newPhi, newThetaDir, newPhiDir, endRad*pow(dr,1.0/2.6), newLen*pow(dr,1.0/2.6), contProb, gam, mass, rgen, rNormGen, vm);
	  newThetaDir = thetaDir + gam*rNormGen->GetValue();
	  rNormGen->Next();
	  if(rgen->GetValue()<0.5){
	    newThetaDir += pi/5.0;
	  } else {
	    newThetaDir -= pi/5.0;
	  }
	  rgen->Next();
	  newPhiDir = phiDir + gam*rNormGen->GetValue();
	  rNormGen->Next();
	  createBranch(newR, newTheta, newPhi, newThetaDir, newPhiDir, endRad*pow(1-dr,1.0/2.6), newLen*pow(1-dr,1.0/2.6), contProb, gam, mass, rgen, rNormGen, vm);
	}
      }
    }
  }
}
			
		
		
				
				
				
							

