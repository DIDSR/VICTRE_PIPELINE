/*
 * duct_options.cxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */


#include "duct.hxx"

using namespace std;
namespace po = boost::program_options;

// default constructor for ductTree
ductTree::ductTree(po::variables_map o, ductTreeInit *init):
  randGen(init->seed),
  u01(randGen),
  radiusDist(o["ductSeg.radiusBetaA"].as<double>(),o["ductSeg.radiusBetaB"].as<double>()){

  opt = o;

  // assign id and update number of ducts
#pragma omp critical
  {
    id = num;
    num += 1;
  }

  fill = vtkImageData::New();
  double spacing[3];
  for(int i=0; i<3; i++){
    spacing[i] = (init->endPos[i] - init->startPos[i])/(init->nFill[i]);
  }
  fill->SetSpacing(spacing);
  fill->SetExtent(0, init->nFill[0]-1, 0, init->nFill[1]-1, 0, init->nFill[2]-1);
  double origin[3];
  for(int i=0; i<3; i++){
    origin[i] = init->startPos[i]+spacing[i]/2.0;
  }
  fill->SetOrigin(origin);

#if VTK_MAJOR_VERSION <= 5
  fill->SetNumberOfScalarComponents(1);
  fill->SetScalarTypeToDouble();
  fill->AllocateScalars();
#else
  fill->AllocateScalars(VTK_DOUBLE,1);
#endif
  
  numBranch = 0;
  maxBranch = o["ductTree.maxBranch"].as<uint>();
  baseLength = o["ductTree.baseLength"].as<double>();

  boundBox = init->boundBox;
  compartmentId = init->compartmentId;
  tissue = init->tissue;
  for(int i=0; i<3; i++){
    prefDir[i] = init->prefDir[i];
  }
  breast = init->breast;
  TDLUloc = init->TDLUloc;
  TDLUattr = init->TDLUattr;

  // temporarily set head branch pointer
  head = nullptr;
}

// destructor
ductTree::~ductTree(){
  fill->Delete();
}

// constructor for first branch (the root)
ductBr::ductBr(double* spos, double* sdir, double r, ductTree *owner){

  double pos[3];
  int invox[3];
  
  bool failSeg = false;	// failed to create a valid segment
  bool edgeSeg = false;	// segment at boundary of ROI

  for(int i=0; i<3; i++){
    startPos[i] = spos[i];
    startDir[i] = sdir[i];
  }
  startRad = r;
  myTree = owner;

  // no parent or sibling branches
  parent = nullptr;
  sibBranch = nullptr;
	
  // root branch has id 0 and level 0 and generation 0
  id = 0;
  level = 0;
  gen = 0;

  // increment tree branch count
  myTree->numBranch += 1;

  // determine length of branch
  length = setLength();
  curLength = 0.0;

  // generate segments to fill branch
  firstSeg = new ductSeg(this);
  lastSeg = firstSeg;

  // update length
  curLength += firstSeg->length;

  if (firstSeg->length == 0.0){
    failSeg = true;
    // debug
    std::cout << "First segment length 0\n" ;
  }

  // check if at ROI boundary by seeing if any neighboring voxels are outside ROI
  double* thePos = lastSeg->endPos;
  double pcoords[3];
  myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
  for(int a=-1; a<=1; a++){
    for(int b=-1; b<=1; b++){
      for(int c=-1; c<=1; c++){
	unsigned char* p =
	  static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	if(p[0] != myTree->compartmentId && p[0] != myTree->tissue->duct){
	  edgeSeg = true;
	  // debug
	  std::cout << "First segment hits boundary\n";
	}
      }
    }
  }

  // generate more segments until proper length
  while(curLength < length && !failSeg && !edgeSeg){
    lastSeg->nextSeg = new ductSeg(lastSeg);  // the lastSeg in parenthesis is used to fill variables including prevSeg ptr
    lastSeg = lastSeg->nextSeg;
    curLength += lastSeg->length;
    if(lastSeg->length == 0.0){
      failSeg = true;
      // debug
      std::cout << "Zero length segment created\n";
    }
    // check if at ROI boundary
    thePos = lastSeg->endPos;
    myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
    for(int a=-1; a<=1; a++){
      for(int b=-1; b<=1; b++){
	for(int c=-1; c<=1; c++){
	  unsigned char* p =
	    static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	  if(p[0] != myTree->compartmentId && p[0] != myTree->tissue->duct){
	    edgeSeg = true;
	    std::cout << "A segment hit the boundary\n";
	  }
	}
      }
    }
  }

  // fill in end of branch variables
  for(int i=0; i<3; i++){
    endPos[i] = lastSeg->endPos[i];
    endDir[i] = lastSeg->endDir[i];
  }
  endRad = lastSeg->endRad;

  length = curLength;

  // set number of children and generate them
  nChild = setChild();

  if(failSeg){
    nChild = 0;
    std::cout << "Segment generation failure for branch" << id << std::endl;
  }

  if(edgeSeg){
    // debug
    //nChild = 0;
    std::cout << "ROI edge collision for branch " << id << std::endl;
  }

  if (nChild == 0){
    firstChild = nullptr;
    secondChild = nullptr;
    // TDLU creation
		
    // check branch length is long enough
    if(length >= myTree->opt["TDLU.minLength"].as<double>()){
      // long enough
      //std::cout << "Adding TDLU" << std::endl;
      
      // pick sizes
      double minLen = myTree->opt["TDLU.minLength"].as<double>();
      double maxLen = myTree->opt["TDLU.maxLength"].as<double>();
      double minWid = myTree->opt["TDLU.minWidth"].as<double>();
      double maxWid = myTree->opt["TDLU.maxWidth"].as<double>();
			
      if(length < maxLen){
	maxLen = length;
      }
			
      double len = minLen + (maxLen-minLen)*myTree->u01();
      double wid = minWid + (maxWid-minWid)*myTree->u01(); 
			
      // save position
      myTree->TDLUloc->InsertNextPoint(endPos);
			
      // save attributes (length width and principle direction)
      double att[5];
      att[0] = len;
      att[1] = wid;
      for(int j=0; j<3; j++){
	att[j+2] = endDir[j];
      }
      myTree->TDLUattr->InsertNextTuple(att);
			
      // segment TDLU
      vtkVector3d axis[3];
      vtkVector3d v2;
      double innerProd;

      // coordinate system
      // first vector
      for(int j=0; j<3; j++){
	axis[0][j] = endDir[j];
      }
			
      // calculate second vector based on direction to coordinate origin
      for(int j=0; j<3; j++){
	v2[j] = endPos[j];
      }
      innerProd = v2.Dot(axis[0]);
      
      for(int j=0; j<3; j++){
	axis[1][j] = v2[j] - innerProd*axis[0][j];
      }
      axis[1].Normalize();

      // calculate 3rd vector based on cross product
      axis[2] = axis[0].Cross(axis[1]);

      // have 3 unit vectors

      double imgRes = myTree->opt["base.imgRes"].as<double>();
      int searchRad = (int)(ceil(len/imgRes));
			
#pragma omp parallel for collapse(3)
      for(int a=-searchRad; a<=searchRad; a++){
	for(int b=-searchRad; b<=searchRad; b++){
	  for(int c=-searchRad; c<=searchRad; c++){
	    double curPos[3] = {endPos[0]+a*imgRes, endPos[1]+b*imgRes, endPos[2]+c*imgRes};
	    int index[3];
	    double pcoords[3];
						
	    // structure coordinates, check if in breast
	    if(myTree->breast->ComputeStructuredCoordinates(curPos, index, pcoords)){
	      // in breast extent
	      
	      // get tissue type, if not duct,skin,nipple,TDLU,outside breast: proceed
	      unsigned char* p = static_cast<unsigned char*>(myTree->breast->GetScalarPointer(index));
							
	      if(p[0] != myTree->tissue->bg && p[0] != myTree->tissue->skin && p[0] != myTree->tissue->nipple && 
		 p[0] != myTree->tissue->TDLU && p[0] != myTree->tissue->duct){
								
		// check if in oval
		// compute position in local coordinate system
		vtkVector3d rvec;
		vtkVector3d lCoords;
		for(int m=0; m<3; m++){
		  rvec[m] = curPos[m]-endPos[m];
		}
		for(int m=0; m<3; m++){
		  lCoords[m] = rvec.Dot(axis[m]);
		}
						
		// inside oval?
		if(lCoords[0]*lCoords[0]/len/len+lCoords[1]*lCoords[1]/wid/wid+lCoords[2]*lCoords[2]/wid/wid < 1.0){
		  p[0] = myTree->tissue->TDLU;
		}
	      }
	    }
	  }
	}
      }
    }
  } else {
    // bifurcate
    // pick radii
    double radii[2];
    double thetas[2];
    setRadiiThetas(radii,thetas);
    // setup first child with level equal to current level
    firstChild = new ductBr(this,level,gen+1,radii[0],thetas[0]);
    firstChild->sibBranch = new ductBr(this,firstChild,level+1,gen+1,radii[1],thetas[1]);
    secondChild = firstChild->sibBranch;
    //for(int i=1; i < nChild; i++){
      // set up other children with higher level than current
    //  lastChild->nextBranch = new ductBr(this,lastChild,level+1,gen+1,radii[i]);
    //  lastChild = lastChild->nextBranch;
    //}
    //free(radii);
  }
}

// constructor for first child branch of a parent branch
ductBr::ductBr(ductBr* par, unsigned int lev, unsigned int g, double r, double theta){

  bool failSeg = false;
  bool edgeSeg = false;
  int invox[3];

  // pointers
  parent = par;
  sibBranch = nullptr;
  //prevBranch = nullptr;
  //nextBranch = this;

  for(int i=0; i<3; i++){
    startPos[i] = parent->endPos[i];
  }
  startRad = r;
  level = lev;
  gen = g;
  myTree = parent->myTree;
  id = myTree->numBranch;
  myTree->numBranch += 1;

  // set starting direction same a parent endDir for first child
  //for(int i=0; i<3; i++){
  //startDir[i] = parent->endDir[i];
  //}
  setDir(startDir,theta);

  // determine length of branch
  length = setLength();
  curLength = 0.0;

  // generate segments to fill branch
  firstSeg = new ductSeg(this);
  lastSeg = firstSeg;

  // update length
  curLength += firstSeg->length;

  if (firstSeg->length == 0.0){
    failSeg = true;
  }

  // check if at ROI boundary
  double* thePos = lastSeg->endPos;
  double pcoords[3];
  myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
  for(int a=-1; a<=1; a++){
    for(int b=-1; b<=1; b++){
      for(int c=-1; c<=1; c++){
	unsigned char* p =
	  static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	if(p[0] != myTree->compartmentId && p[0] != myTree->tissue->duct){
	  edgeSeg = true;
	}
      }
    }
  }

  // generate more segments until proper length
  while(curLength < length && !failSeg && !edgeSeg){
    lastSeg->nextSeg = new ductSeg(lastSeg);  // the lastSeg in parenthesis is used to fill variables including prevSeg ptr
    lastSeg = lastSeg->nextSeg;
    curLength += lastSeg->length;
    if(lastSeg->length == 0.0){
      failSeg = true;
    }
    // check if at ROI boundary
    thePos = lastSeg->endPos;
    myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
    for(int a=-1; a<=1; a++){
      for(int b=-1; b<=1; b++){
	for(int c=-1; c<=1; c++){
	  unsigned char* p =
	    static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	  if(p[0] != myTree->compartmentId && p[0] != myTree->tissue->duct){
	    edgeSeg = true;
	  }
	}
      }
    }
  }

  // fill in end of branch variables
  for(int i=0; i<3; i++){
    endPos[i] = lastSeg->endPos[i];
    endDir[i] = lastSeg->endDir[i];
  }
  endRad = lastSeg->endRad;

  length = curLength;

  // set number of children and generate them
  nChild = setChild();

  if(failSeg){
    nChild = 0;
    std::cout << "Segment generation failure for branch" << id << std::endl;
  }

  if(edgeSeg){
    nChild = 0;
    //std::cout << "ROI edge collision for branch " << id << std::endl;
  }

  if (nChild == 0){
    firstChild = nullptr;
    secondChild = nullptr;
		
    // TDLU creation
		
    // check branch length is long enough
    if(length >= myTree->opt["TDLU.minLength"].as<double>()){
      // long enough
      //std::cout << "Adding TDLU" << std::endl;
		
      // pick sizes
      double minLen = myTree->opt["TDLU.minLength"].as<double>();
      double maxLen = myTree->opt["TDLU.maxLength"].as<double>();
      double minWid = myTree->opt["TDLU.minWidth"].as<double>();
      double maxWid = myTree->opt["TDLU.maxWidth"].as<double>();
			
      if(length < maxLen){
	maxLen = length;
      }
			
      double len = minLen + (maxLen-minLen)*myTree->u01();
      double wid = minWid + (maxWid-minWid)*myTree->u01(); 
			
      // save position
      myTree->TDLUloc->InsertNextPoint(endPos);
			
      // save attributes (length width and principle direction)
      double att[5];
      att[0] = len;
      att[1] = wid;
      for(int j=0; j<3; j++){
	att[j+2] = endDir[j];
      }
      myTree->TDLUattr->InsertNextTuple(att);
			
      // segment TDLU
      vtkVector3d axis[3];
      vtkVector3d v2;
      double innerProd;

      // coordinate system
      // first vector
      for(int j=0; j<3; j++){
	axis[0][j] = endDir[j];
      }
			
      // calculate second vector based on direction to coordinate origin
      for(int j=0; j<3; j++){
	v2[j] = endPos[j];
      }
      innerProd = v2.Dot(axis[0]);

      for(int j=0; j<3; j++){
	axis[1][j] = v2[j] - innerProd*axis[0][j];
      }
      axis[1].Normalize();

      // calculate 3rd vector based on cross product
      axis[2] = axis[0].Cross(axis[1]);

      // have 3 unit vectors

      double imgRes = myTree->opt["base.imgRes"].as<double>();
      int searchRad = (int)(ceil(len/imgRes));
			
#pragma omp parallel for collapse(3)
      for(int a=-searchRad; a<=searchRad; a++){
	for(int b=-searchRad; b<=searchRad; b++){
	  for(int c=-searchRad; c<=searchRad; c++){
	    double curPos[3] = {endPos[0]+a*imgRes, endPos[1]+b*imgRes, endPos[2]+c*imgRes};
	    int index[3];
	    double pcoords[3];
						
	    // structure coordinates, check if in breast
	    if(myTree->breast->ComputeStructuredCoordinates(curPos, index, pcoords)){
	      // in breast extent
	      
	      // get tissue type, if not duct,skin,nipple,TDLU,outside breast: proceed
	      unsigned char* p = static_cast<unsigned char*>(myTree->breast->GetScalarPointer(index));
							
	      if(p[0] != myTree->tissue->bg && p[0] != myTree->tissue->skin && p[0] != myTree->tissue->nipple && 
		 p[0] != myTree->tissue->TDLU && p[0] != myTree->tissue->duct){
								
		// check if in oval
		// compute position in local coordinate system
		vtkVector3d rvec;
		vtkVector3d lCoords;
		for(int m=0; m<3; m++){
		  rvec[m] = curPos[m]-endPos[m];
		}
		for(int m=0; m<3; m++){
		  lCoords[m] = rvec.Dot(axis[m]);
		}
						
		// inside oval?
		if(lCoords[0]*lCoords[0]/len/len+lCoords[1]*lCoords[1]/wid/wid+lCoords[2]*lCoords[2]/wid/wid < 1.0){
		  p[0] = myTree->tissue->TDLU;
		}
	      }
	    }
	  }
	}
      }
    }
  } else {
    // pick radii and thetas
    double radii[2];
    double thetas[2];
    setRadiiThetas(radii,thetas);
    // setup first child with level equal to current level
    firstChild = new ductBr(this,level,gen+1,radii[0],thetas[0]);
    firstChild->sibBranch = new ductBr(this,firstChild,level+1,gen+1,radii[1],thetas[1]);
    secondChild = firstChild->sibBranch;    
    //lastChild = firstChild;
    //for(int i=1; i < nChild; i++){
      // set up other children with higher level than current
    //  lastChild->nextBranch = new ductBr(this,lastChild,level+1,gen+1,radii[i]);
    //  lastChild = lastChild->nextBranch;
    //}
    //free(radii);
  }
}

// constructor for subsequent children (not first child) of a parent branch
ductBr::ductBr(ductBr* par, ductBr* par2, unsigned int lev, unsigned int g, double r, double theta){

  bool failSeg = false;
  bool edgeSeg = false;

  int invox[3];

  // pointers
  parent = par;
  sibBranch = par2;
  //nextBranch = this;

  for(int i=0; i<3; i++){
    startPos[i] = parent->endPos[i];
  }
  startRad = r;
  level = lev;
  gen = g;
  myTree = parent->myTree;
  id = myTree->numBranch;
  myTree->numBranch += 1;

  // set starting direction
  setDir(startDir,theta);

  // determine length of branch
  length = setLength();
  curLength = 0.0;

  // generate segments to fill branch
  firstSeg = new ductSeg(this);
  lastSeg = firstSeg;

  // update length
  curLength += firstSeg->length;

  if (firstSeg->length == 0.0){
    failSeg = true;
  }

  // check if at ROI boundary
  double* thePos = lastSeg->endPos;
  double pcoords[3];
  myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
  for(int a=-1; a<=1; a++){
    for(int b=-1; b<=1; b++){
      for(int c=-1; c<=1; c++){
	unsigned char* p =
	  static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	if(p[0] != myTree->compartmentId && p[0] != myTree->tissue->duct){
	  edgeSeg = true;
	}
      }
    }
  }

  // generate more segments until proper length
  while(curLength < length && !failSeg && !edgeSeg){
    lastSeg->nextSeg = new ductSeg(lastSeg);  // the lastSeg in parenthesis is used to fill variables including prevSeg ptr
    lastSeg = lastSeg->nextSeg;
    curLength += lastSeg->length;
    if(lastSeg->length == 0.0){
      failSeg = true;
    }
    // check if at ROI boundary
    thePos = lastSeg->endPos;
    myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
    for(int a=-1; a<=1; a++){
      for(int b=-1; b<=1; b++){
	for(int c=-1; c<=1; c++){
	  unsigned char* p =
	    static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	  if(p[0] != myTree->compartmentId && p[0] != myTree->tissue->duct){
	    edgeSeg = true;
	  }
	}
      }
    }
  }

  // fill in end of branch variables
  for(int i=0; i<3; i++){
    endPos[i] = lastSeg->endPos[i];
    endDir[i] = lastSeg->endDir[i];
  }
  endRad = lastSeg->endRad;

  length = curLength;

  // set number of children and generate them
  nChild = setChild();

  if(failSeg){
    nChild = 0;
    std::cout << "Segment generation failure for branch" << id << std::endl;
  }

  if(edgeSeg){
    nChild = 0;
    //std::cout << "ROI edge collision for branch " << id << std::endl;
  }

  if (nChild == 0){
    firstChild = nullptr;
    secondChild = nullptr;
		
    // TDLU creation
		
    // check branch length is long enough
    if(length >= myTree->opt["TDLU.minLength"].as<double>()){
      // long enough
      //std::cout << "Adding TDLU" << std::endl;
		
      // pick sizes
      double minLen = myTree->opt["TDLU.minLength"].as<double>();
      double maxLen = myTree->opt["TDLU.maxLength"].as<double>();
      double minWid = myTree->opt["TDLU.minWidth"].as<double>();
      double maxWid = myTree->opt["TDLU.maxWidth"].as<double>();
			
      if(length < maxLen){
	maxLen = length;
      }
			
      double len = minLen + (maxLen-minLen)*myTree->u01();
      double wid = minWid + (maxWid-minWid)*myTree->u01(); 
			
      // save position
      myTree->TDLUloc->InsertNextPoint(endPos);
			
      // save attributes (length width and principle direction)
      double att[5];
      att[0] = len;
      att[1] = wid;
      for(int j=0; j<3; j++){
	att[j+2] = endDir[j];
      }
      myTree->TDLUattr->InsertNextTuple(att);
      
      // segment TDLU
      vtkVector3d axis[3];
      vtkVector3d v2;
      double innerProd;

      // coordinate system
      // first vector
      for(int j=0; j<3; j++){
	axis[0][j] = endDir[j];
      }
			
      // calculate second vector based on direction to coordinate origin
      for(int j=0; j<3; j++){
	v2[j] = endPos[j];
      }
      innerProd = v2.Dot(axis[0]);
      
      for(int j=0; j<3; j++){
	axis[1][j] = v2[j] - innerProd*axis[0][j];
      }
      axis[1].Normalize();
      
      // calculate 3rd vector based on cross product
      axis[2] = axis[0].Cross(axis[1]);

      // have 3 unit vectors

      double imgRes = myTree->opt["base.imgRes"].as<double>();
      int searchRad = (int)(ceil(len/imgRes));
			
#pragma omp parallel for collapse(3)
      for(int a=-searchRad; a<=searchRad; a++){
	for(int b=-searchRad; b<=searchRad; b++){
	  for(int c=-searchRad; c<=searchRad; c++){
	    double curPos[3] = {endPos[0]+a*imgRes, endPos[1]+b*imgRes, endPos[2]+c*imgRes};
	    int index[3];
	    double pcoords[3];
						
	    // structure coordinates, check if in breast
	    if(myTree->breast->ComputeStructuredCoordinates(curPos, index, pcoords)){
	      // in breast extent
	      
	      // get tissue type, if not duct,skin,nipple,TDLU,outside breast: proceed
	      unsigned char* p = static_cast<unsigned char*>(myTree->breast->GetScalarPointer(index));
							
	      if(p[0] != myTree->tissue->bg && p[0] != myTree->tissue->skin && p[0] != myTree->tissue->nipple && 
		 p[0] != myTree->tissue->TDLU && p[0] != myTree->tissue->duct){
		
		// check if in oval
		// compute position in local coordinate system
		vtkVector3d rvec;
		vtkVector3d lCoords;
		for(int m=0; m<3; m++){
		  rvec[m] = curPos[m]-endPos[m];
		}
		for(int m=0; m<3; m++){
		  lCoords[m] = rvec.Dot(axis[m]);
		}
						
		// inside oval?
		if(lCoords[0]*lCoords[0]/len/len+lCoords[1]*lCoords[1]/wid/wid+lCoords[2]*lCoords[2]/wid/wid < 1.0){
		  p[0] = myTree->tissue->TDLU;
		}
	      }
	    }
	  }
	}
      }
    }
		
  } else {
    // pick radii
    double radii[2];
    double thetas[2];
    setRadiiThetas(radii,thetas);

    // setup first child with level equal to current level
    firstChild = new ductBr(this,level,gen+1,radii[0],thetas[0]);
    firstChild->sibBranch = new ductBr(this,firstChild,level+1,gen+1,radii[1],thetas[1]);
    secondChild = firstChild->sibBranch;
    //lastChild = firstChild;
    //for(int i=1; i < nChild; i++){
      // set up other children with higher level than current
    //  lastChild->nextBranch = new ductBr(this,lastChild,level+1,gen+1,radii[i]);
    //  lastChild = lastChild->nextBranch;
    //}
    //free(radii);
  }
}

double ductBr::setLength(void){
  // set length using random distribution and level
  double len;
  double randVal = myTree->u01();
  double baseLen = myTree->baseLength;
  double lenShrink = myTree->opt["ductBr.lenShrink"].as<double>();
  double lenRange = myTree->opt["ductBr.lenRange"].as<double>();

  len = baseLen*pow(lenShrink,level);

  // add variability +/- lenRange fraction of len
  len = len - lenRange*len + randVal*2*lenRange*len;
  return(len);
}

unsigned int ductBr::setChild(void){
  // determine number of child branches
  // cumulative probabilities of having 0-4 children based on level
  //unsigned int maxChild = myTree->opt["ductBr.maxChild"].as<uint>();
  //unsigned int levBound = myTree->opt["ductBr.childLevBound"].as<uint>();
  //ostringstream var;

  //double *prob = (double *)malloc((levBound+1)*maxChild*sizeof(double));

  //for(int a=0; a<=levBound; a++){
  //  for(int b=0; b<maxChild; b++){
  //    var.str("");
  //    var.clear();
  //    var << "ductBr.child" << a << b;
  //    prob[a*maxChild+b] = myTree->opt[var.str()].as<double>();
  //  }
  //}

  // if small enough, no children
  double minRad = myTree->opt["ductBr.childMinRad"].as<double>();
  if(endRad < minRad){
    //free(prob);
    return(0);
  }

  // check for max number branches
  if(myTree->numBranch >= myTree->maxBranch){
    //free(prob);
    return(0);
  }

  // define maximum generation
  unsigned int maxGen = myTree->opt["ductTree.maxGen"].as<uint>();
  if(gen > maxGen){
    //free(prob);
    return(0);
  }

  //double randVal = myTree->u01();
  //
  //for(int b=maxChild-1; b >= 0; b--){
  //  if(randVal > prob[std::min(level,levBound)*maxChild+b]){
  //    free(prob);
  //    return(b+1);
  //  }
  //}
  // default to no children
  //free(prob);

  // default 2 children
  return(2);
}

void ductBr::setRadiiThetas(double* radii, double* thetas){
  // set radii and angles of child branches based on the parent
  
  double minFrac = myTree->opt["ductBr.minRadFrac"].as<double>();
  double maxFrac = myTree->opt["ductBr.maxRadFrac"].as<double>();

  double randVal = myTree->u01();

  // first child radius
  double myFrac = minFrac + randVal*(maxFrac-minFrac);
  radii[0] = myFrac*endRad;

  // second child radius based on Murray's law
  double b = 1.0/myFrac;
  double a = pow(pow(b,3.0)-1.0,1.0/3.0);
  radii[1] = a*radii[0]; 

  // first child theta
  double lBound = (pow(b,4.0)+1.0-pow(a,4.0))/(2.0*pow(b,2.0));
  double uBound = (pow(b,2.0)+1.0-pow(a,2.0))/(2.0*b);
  randVal = myTree->u01();
  double ctheta = lBound + randVal*(uBound-lBound);
  thetas[0] = acos(ctheta);

  // second child theta
  lBound = (pow(b,4.0)+pow(a,4.0)-1.0)/(2*pow(a,2.0)*pow(b,2.0));
  uBound = (pow(b,2.0)+pow(a,2.0)-1.0)/(2*a*b);
  randVal = myTree->u01();
  ctheta = lBound + randVal*(uBound-lBound);
  thetas[1] = acos(ctheta);
}

void ductBr::setDir(double* sdir, double theta){

  const double pi = boost::math::constants::pi<double>();

  // set initial direction of branch
  double dir[3];
  double tempV[3];
  double basis1[3];
  double basis2[3];

  double rotateJitter = myTree->opt["ductBr.rotateJitter"].as<double>(); 
  double rotate;

  if(sibBranch == nullptr){
    // this is the first child
    // random rotation about parent direction
    rotate = 2*pi*myTree->u01();
    azimuth = rotate;
  } else {
    // this is the second child
    rotate = sibBranch->azimuth + pi;
    double randVal = myTree->u01();
    rotate = rotate - rotateJitter + randVal*2*rotateJitter;
    azimuth = rotate;
  } 

  // project origin onto plane perpendicular parent endDir
  double dotProd = 0.0;
  for(int i=0; i<3; i++){
    dotProd += parent->endDir[i]*startPos[i];
  }
  for(int i=0; i<3; i++){
    tempV[i] = dotProd*parent->endDir[i];
    basis1[i] = tempV[i] - startPos[i];
  }

  // normalize basis1
  double norm = 0.0;
  for(int i=0; i<3; i++){
    norm += basis1[i]*basis1[i];
  }
  norm = sqrt(norm);
  for(int i=0; i<3; i++){
    basis1[i] = basis1[i]/norm;
  }

  // find second basis vector using cross product
  basis2[0] = parent->endDir[1]*basis1[2] - parent->endDir[2]*basis1[1];
  basis2[1] = parent->endDir[2]*basis1[0] - parent->endDir[0]*basis1[2];
  basis2[2] = parent->endDir[0]*basis1[1] - parent->endDir[1]*basis1[0];

  for(int i=0; i<3; i++){
    sdir[i] = cos(theta)*parent->endDir[i] + sin(theta)*(cos(rotate)*basis1[i] + sin(rotate)*basis2[i]);
  }
}

// branch destructor that also deletes all child branches
ductBr::~ductBr(){
  // delete segments
  ductSeg* delSeg;
  while(firstSeg != lastSeg){
    delSeg = firstSeg;
    firstSeg = firstSeg->nextSeg;
    delete(delSeg);
  }
  delete(firstSeg);

  // delete child branches
  ductBr* delBranch;
  while(firstChild != secondChild){
    delBranch = firstChild;
    firstChild = firstChild->secondChild;
    delete(delBranch);
  }
  delete(firstChild);
}

// constructor for first segment
ductSeg::ductSeg(ductBr* br){
  myBranch = br;
  prevSeg = nullptr;
  nextSeg = this;

  for(int i=0; i<3; i++){
    startPos[i] = myBranch->startPos[i];
    startDir[i] = myBranch->startDir[i];
  }
  startRad = myBranch->startRad;

  // keeping derivatives zero at nodes for now
  startDeriv = 0.0;

  // code to generate random segment
  makeSeg();
}

// constructor for subsequent segments
ductSeg::ductSeg(ductSeg* pr){
  prevSeg = pr;
  myBranch = prevSeg->myBranch;
  nextSeg = this;

  for(int i=0; i<3; i++){
    startPos[i] = prevSeg->endPos[i];
    startDir[i] = prevSeg->endDir[i];
  }

  startRad = prevSeg->endRad;
  startDeriv = prevSeg->endDeriv;

  // code to generate random segment
  makeSeg();
}

void ductSeg::makeSeg(){

  const double pi = boost::math::constants::pi<double>();

  double segFrac = myBranch->myTree->opt["ductSeg.segFrac"].as<double>();
  //double minLen = myBranch->myTree->opt["ductSeg.minLen"].as<double>();
  //double maxLen = myBranch->myTree->opt["ductSeg.maxLen"].as<double>();
  unsigned int numTry = myBranch->myTree->opt["ductSeg.numTry"].as<uint>();
  unsigned int maxTry = myBranch->myTree->opt["ductSeg.maxTry"].as<uint>();
  unsigned int absMaxTry = myBranch->myTree->opt["ductSeg.absMaxTry"].as<uint>();
  double maxRad = myBranch->myTree->opt["ductSeg.maxCurvRad"].as<double>();
  double angleMax =  pi*myBranch->myTree->opt["ductSeg.maxCurvFrac"].as<double>();
  double roiStep = myBranch->myTree->opt["ductSeg.roiStep"].as<double>();
  double densityWt = myBranch->myTree->opt["ductSeg.densityWt"].as<double>();
  double angleWt = myBranch->myTree->opt["ductSeg.angleWt"].as<double>();
  double prefDir[3]; // preferential direction of growth
  for(int i=0; i<3; i++){
    prefDir[i] = myBranch->myTree->prefDir[i];
  }
  double maxEndRad = myBranch->myTree->opt["ductSeg.maxEndRad"].as<double>();
  double minEndRad = myBranch->myTree->opt["ductSeg.minEndRad"].as<double>();

  int fillExtent[6];	// extents of fill
  myBranch->myTree->fill->GetExtent(fillExtent);

  double pos[3];
  unsigned int invox[3];

  unsigned int curTry;		// number of valid segments tested so far
  unsigned int totalTry;		// number of test segments so far
  unsigned int allTry;		// total number of test segments overall

  double lengthLB, lengthUB;  // bounds on random length
  double randVal, quantileVal;  // for length random generator

  double theta;  				// rotation of segment
  double radius; 				// segment radius of curvature
  double radLB,radUB; 		// min and max radius of curvature

  double curv[3];  			// point of rotation
  double curvNorm;			// stores norm(startPos-curv)
  double basis1[3];
  double basis2[3]; 			// basis vectors
  double tempV[3];
  double checkPos[3];			// checking if segment position in ROI

  int myVoxel[3];
  double pcoords[3];
  double checkAngle;
  double checkLength;
  double currentDir[3];  // unit vector in current direction
  double angleStep; // angular step size for checking in ROI for segment

  bool foundSeg = false;  // have we found a good segment yet?
  bool inROI;		// is current test segment in ROI?
  bool inFOV;		// is current segment in FOV?

  // breast FOV for checking if segment in FOV
  double breastOrigin[3];
  double breastSpacing[3];
  int breastDim[3];
  double breastFOV[6];

  myBranch->myTree->breast->GetOrigin(breastOrigin);
  myBranch->myTree->breast->GetSpacing(breastSpacing);
  myBranch->myTree->breast->GetDimensions(breastDim);

  for(int i=0; i<3; i++){
    breastFOV[2*i] = breastOrigin[i];
    breastFOV[2*i+1] = breastOrigin[i]+(double)breastDim[i]*breastSpacing[i];
  }

  double bestCurv[3];	// best center of curvature found so far
  double bestRadius;
  double bestCost;  // best cost found so far
  double cost;		// current cost

  double density;	// change in density of duct structure due to new segment

  // determine proposed segment length
  // default length
  length = myBranch->length*segFrac;

  if(length > (myBranch->length - myBranch->curLength)){
    // truncate
    length = myBranch->length - myBranch->curLength;
  }
  //if (maxLen < (myBranch->length - myBranch->curLength)/10.0){
  //  length = maxLen;
  //} else if (minLen > (myBranch->length - myBranch->curLength)) {
  //  length = minLen;
  //} else {
    // use beta distribution to choose random length
    // beta distribution bounds;
  //  lengthLB = fmax(minLen,(myBranch->length - myBranch->curLength)/10.0);
  //  lengthUB = fmin(maxLen,(myBranch->length - myBranch->curLength));

    // generate a random length using beta distribution
  //  randVal = myBranch->myTree->u01();
  //  quantileVal = boost::math::quantile(myBranch->myTree->lengthDist, randVal);
    // scale to length range
  //  length = quantileVal*(lengthUB-lengthLB) + lengthLB;
  //}

  allTry = 0;

  while (!foundSeg && allTry < absMaxTry){
    curTry = 0;
    while (curTry < numTry){
      totalTry = 0;
      inROI = false;
      inFOV = false;
      while (!inROI && !inFOV && totalTry < maxTry){
	allTry++;
	// generate random segment
	theta = 2*pi*myBranch->myTree->u01();
	radUB = maxRad;
	radLB = length/angleMax;
	// use beta distribution to pick radius
	randVal = myBranch->myTree->u01();
	quantileVal = boost::math::quantile(myBranch->myTree->radiusDist, randVal);
	// scale to radius range
	radius = quantileVal*(radUB-radLB) + radLB;
	totalTry += 1;

	// checking if in ROI
	// need basis vectors in plane perpendicular to startDir

	// project origin (0,0,0) onto plane perpendicular to startDir
	vtkMath::ProjectVector(startPos, startDir, tempV);
	
	vtkMath::Subtract(tempV, startPos, basis1);

	// normalize it
	vtkMath::Normalize(basis1);

	// find second basis vector using cross product
	vtkMath::Cross(startDir,basis1,basis2);

	// calculate curvature
	for(int i=0; i<3; i++){
	  curv[i] = startPos[i] + radius*(basis1[i]*cos(theta) + basis2[i]*sin(theta));
	}

	// calculate norm(startPos-curv)
	curvNorm = 0.0;
	for(int i=0; i<3; i++){
	  curvNorm += (startPos[i]-curv[i])*(startPos[i]-curv[i]);
	}
	curvNorm = sqrt(curvNorm);

	// check if in ROI
	angleStep = roiStep/radius;
	checkAngle = 0.0;
	checkLength = 0.0;
	inROI = true;
	inFOV = true;
	while (checkLength < length && inROI && inFOV){
	  for(int i=0; i<3; i++){
	    checkPos[i] = curv[i] + radius*((startPos[i]-curv[i])/curvNorm*cos(checkAngle)+startDir[i]*sin(checkAngle));
	  }

	  // is point in FOV and in ROI?
	  
	  // check FOV first
	  if(checkPos[0] < breastFOV[0] || checkPos[0] > breastFOV[1] ||
	     checkPos[1] < breastFOV[2] || checkPos[1] > breastFOV[3] ||
	     checkPos[2] < breastFOV[4] || checkPos[2] > breastFOV[5]){
	    inFOV = false;
	  }

	  // check in ROI
	  
	  if(inFOV){
	    myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, myVoxel, pcoords);

	    unsigned char* p =
	      static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(myVoxel[0],myVoxel[1],myVoxel[2]));
	    if(p[0] != myBranch->myTree->compartmentId && p[0] != myBranch->myTree->tissue->duct){
	      inROI = false;
	    }
	  }

	  checkAngle += angleStep;
	  checkLength += angleStep*radius;
	}
	// check the end point
	for(int i=0; i<3; i++){
	  checkPos[i] = curv[i] + radius*((startPos[i]-curv[i])/curvNorm*cos(length/radius)+startDir[i]*sin(length/radius));
	}
	// check FOV first
	if(checkPos[0] < breastFOV[0] || checkPos[0] > breastFOV[1] ||
	   checkPos[1] < breastFOV[2] || checkPos[1] > breastFOV[3] ||
	   checkPos[2] < breastFOV[4] || checkPos[2] > breastFOV[5]){
	  inFOV = false;
	}

	// check in ROI

	if(inFOV){
	  myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, myVoxel, pcoords);

	  unsigned char* p =
	    static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(myVoxel[0],myVoxel[1],myVoxel[2]));
	  if(p[0] != myBranch->myTree->compartmentId && p[0] != myBranch->myTree->tissue->duct){
	    inROI = false;
	  }
	}
      }
      curTry += 1;
      // if valid
      if (inROI && inFOV){
	if(!foundSeg){
	  foundSeg = true;
	  // this is first valid segment, must be the best
	  // calculate cost and set to current best
	  
	  // reduction in squared distance to ducts in ROI
	  // only evaluate endPos
	  density = 0.0;

	  // iterate over fill voxels
#pragma omp parallel for collapse(3) reduction(+:density)
	  for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
	    for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
	      for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
		double* v = static_cast<double*>(myBranch->myTree->fill->GetScalarPointer(a,b,c));
		if(v[0] > 0.0){
		  double dist;
		  // voxel in ROI, calculate change in distance
		  // voxel location
		  vtkIdType id;
		  int coord[3];
		  coord[0] = a;
		  coord[1] = b;
		  coord[2] = c;
		  id = myBranch->myTree->fill->ComputePointId(coord);
		  // get spatial coordinates of point
		  double pos[3];
		  myBranch->myTree->fill->GetPoint(id,pos);
		  dist = vtkMath::Distance2BetweenPoints(checkPos, pos);
		  if(dist < v[0]){
		    density -= v[0] - dist;
		  }
		}
	      }
	    }
	  }

	  // penalty includes direction of segment (away from preferential direction)
	  // negative cost is good, dot product gives cosine of angle
	  // endDir from derivative of position
	  for(int i=0; i<3; i++){
	    endDir[i] = -1*(startPos[i]-curv[i])/curvNorm*sin(length/radius)+startDir[i]*cos(length/radius);
	  }
	  // normalize
	  vtkMath::Normalize(endDir);
	  
	  bestCost = densityWt*density - angleWt*vtkMath::Dot(endDir,prefDir);
	  bestRadius = radius;
	  for(int i=0; i<3; i++){
	    bestCurv[i] = curv[i];
	  }
	} else {
	  // calculate segment cost, if best yet, update current best
	  density = 0.0;

	  // iterate over fill voxels
#pragma omp parallel for collapse(3) reduction(+:density)
	  for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
	    for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
	      for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
		double* v = static_cast<double*>(myBranch->myTree->fill->GetScalarPointer(a,b,c));
		if(v[0] > 0.0){
		  double dist;
		  // voxel in ROI, calculate change in distance
		  // voxel location
		  vtkIdType id;
		  int coord[3];
		  coord[0] = a;
		  coord[1] = b;
		  coord[2] = c;
		  id = myBranch->myTree->fill->ComputePointId(coord);
		  // get spatial coordinates of point
		  double pos[3];
		  myBranch->myTree->fill->GetPoint(id,pos);
		  dist = vtkMath::Distance2BetweenPoints(checkPos, pos);
		  if(dist < v[0]){
		    density -= v[0] - dist;
		  }
		}
	      }
	    }
	  }

	  // endDir from derivative of position
	  for(int i=0; i<3; i++){
	    endDir[i] = -1*(startPos[i]-curv[i])/curvNorm*sin(length/radius)+startDir[i]*cos(length/radius);
	  }
	  // normalize
	  vtkMath::Normalize(endDir);

	  cost = densityWt*density - angleWt*vtkMath::Dot(endDir,prefDir);
	  if (cost < bestCost){
	    // found a new best segment
	    bestCost = cost;
	    bestRadius = radius;
	    for(int i=0; i<3; i++){
	      bestCurv[i] = curv[i];
	    }
	  }
	}
      }
    }
    if(!foundSeg){
      // couldn't find good segment, reduce length
      length = length/10.0;  // could get into infinite loop
    }
  }
  if(!foundSeg){
    // we have failed completely
    length = 0.0;
    for(int i=0; i<3; i++){
      endPos[i] = startPos[i];
      endDir[i] = startDir[i];
    }
    endRad = startRad;
  } else {
    // determined new segment, fill in variables
    curvNorm = 0.0;
    for(int i=0; i<3; i++){
      curvNorm += (startPos[i]-bestCurv[i])*(startPos[i]-bestCurv[i]);
    }
    curvNorm = sqrt(curvNorm);

    for(int i=0; i<3; i++){
      centerCurv[i] = bestCurv[i];
      endPos[i] = centerCurv[i] + bestRadius*((startPos[i]-bestCurv[i])/curvNorm*cos(length/bestRadius)+startDir[i]*sin(length/bestRadius));
    }
    radCurv = 0.0;
    for(int i=0; i<3; i++){
      radCurv += (centerCurv[i]-startPos[i])*(centerCurv[i]-startPos[i]);
    }
    radCurv = sqrt(radCurv);

    // length already set
    // endDir from derivative of position
    for(int i=0; i<3; i++){
      endDir[i] = -1*(startPos[i]-centerCurv[i])/radCurv*sin(length/bestRadius)+startDir[i]*cos(length/bestRadius);
    }
    // normalize
    vtkMath::Normalize(endDir);

    // keeping derivatives fixed to 0 for now
    endDeriv = 0.0;
    // end radius from uniform random variable
    randVal = myBranch->myTree->u01();
    endRad = minEndRad*startRad + randVal*(maxEndRad-minEndRad)*startRad;
    // shape parameters
    setShape();
    // update voxel-based visualization
    updateMap();
    // update fill
#pragma omp parallel for collapse(3)
    for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
      for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
	for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
	  double* v = static_cast<double*>(myBranch->myTree->fill->GetScalarPointer(a,b,c));
	  if(v[0] > 0.0){
	    double dist;
	    // voxel in ROI
	    // voxel location
	    vtkIdType id;
	    int coord[3];
	    coord[0] = a;
	    coord[1] = b;
	    coord[2] = c;
	    id = myBranch->myTree->fill->ComputePointId(coord);
	    // get spatial coordinates of point
	    double pos[3];
	    myBranch->myTree->fill->GetPoint(id,pos);
	    dist = vtkMath::Distance2BetweenPoints(endPos, pos);
	    if(dist < v[0]){
	      // update minimum distance
	      v[0] = dist;
	    }
	  }
	}
      }
    }
  }
}

void ductSeg::setShape(){
  // cubic polynomial v(0)d^3+v(1)d^2+v(2)d+v(3)
  // d = [0,length)
  shape[3] = startRad;
  shape[2] = startDeriv;
  double c = endRad-shape[2]*length-shape[3];
  shape[0] = (endDeriv-shape[2]-2*c/length)/length/length;
  shape[1] = c/length/length-shape[0]*length;
}

double ductSeg::getRadius(double t){
  // return segment radius at position t (mm)
  if(t>=0 && t<=length){
    return(shape[0]*t*t*t+shape[1]*t*t+shape[2]*t+shape[3]);
  } else {
    return(0);
  }
}

void ductSeg::updateMap(){

  const double pi = boost::math::constants::pi<double>();

  // update voxelized map of ducts
  double basis1[3];
  double basis2[3];
  double basis3[3];

  for(int i=0; i<3; i++){
    basis1[i] = startDir[i];
    basis2[i] = (centerCurv[i] - startPos[i])/radCurv;
  }

  vtkMath::Cross(basis1,basis2,basis3);

  // step size for updating, half of voxel width
  double spacing[3];
  myBranch->myTree->breast->GetSpacing(spacing);
  double step = spacing[0];
  if(spacing[1] < step){
    step = spacing[1];
  }
  if(spacing[2] < step){
    step = spacing[2];
  }
  step = step/2.0;
  
  // step along segment, radius
  double ls,rs;
	
  // step along length of segment
  ls = step;
  // step along radius of segment
  rs = step;

  // position along segment
  //double lpos=0.0;

  // calculate number of for loops for openMP
  int lIter = (int)(ceil(length/ls))+1;

#pragma omp parallel for
  for(int j=0; j<=lIter; j++){
    double lpos = j*ls;
    double currentRad = getRadius(lpos);
    double rpos = 0.0;

    double currentPos[3];
    double currentDir[3];
    double lbasis2[3];
		
    // current position
    for(int i=0; i<3; i++){
      currentPos[i] = centerCurv[i] + radCurv*(-1*cos(lpos/radCurv)*basis2[i] + sin(lpos/radCurv)*basis1[i]);
      currentDir[i] = sin(lpos/radCurv)*basis2[i] + cos(lpos/radCurv)*basis1[i];
      lbasis2[i] = (centerCurv[i] - currentPos[i])/radCurv;
    }

    // plane perpendicular to currentDir spanned by basis3 and lbasis2
    while(rpos < currentRad){
      double checkPos[3];
      int checkIdx[3];
      double pcoords[3];
      
      if(rpos < step){
	// only check current voxel assume angle = 0
	for(int i=0; i<3; i++){
	  checkPos[i] = currentPos[i] + rpos*(-1*cos(0.0)*lbasis2[i] + sin(0.0)*basis3[i]);
	}
	myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, checkIdx, pcoords);
	// set voxel to duct
	unsigned char* p =
	  static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(checkIdx[0],checkIdx[1],checkIdx[2]));
	p[0] = myBranch->myTree->tissue->duct;
      } else {
	// angle step
	double as = step/rpos;
	double apos = 0.0;
	while(apos < 2*pi){
	  for(int i=0; i<3; i++){
	    checkPos[i] = currentPos[i] + rpos*(-1*cos(apos)*lbasis2[i] + sin(apos)*basis3[i]);
	  }
	  myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, checkIdx, pcoords);
	  // set voxel to duct
	  unsigned char* p =
	    static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(checkIdx[0],checkIdx[1],checkIdx[2]));
	  p[0] = myBranch->myTree->tissue->duct;
	  apos += as;
	}
      }
      rpos += rs;
    }
    
    //lpos += ls;
  }
}
