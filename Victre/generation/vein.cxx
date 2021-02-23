/*
 * vein_options.cxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#include "vein.hxx"

using namespace std;
namespace po = boost::program_options;

// default constructor for veinTree
veinTree::veinTree(po::variables_map o, veinTreeInit *init):
  randGen(init->seed),
  u01(randGen),
    radiusDist(o["vesselSeg.radiusBetaA"].as<double>(),o["vesselSeg.radiusBetaB"].as<double>()){

  opt = o;

  // assign id and update number of arteries
  id = num;
  num += 1;

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
  maxBranch = o["vesselTree.maxBranch"].as<uint>();
  baseLength = o["vesselTree.baseLength"].as<double>();

  boundBox = init->boundBox;
  tissue = init->tissue;
  for(int i=0; i<3; i++){
    nipplePos[i] = init->nipplePos[i];
  }
  breast = init->breast;

  // temporarily set head branch pointer
  head = nullptr;
}

// destructor
veinTree::~veinTree(){
  fill->Delete();
}


// constructor for first branch (the root)
veinBr::veinBr(double* spos, double* sdir, double r, veinTree *owner){

  double pos[3];
  int invox[3];

  bool failSeg = false;	// failed to create a valid segment
  bool edgeSeg = false;	// segment at boundary of ROI

  myTree = owner;

  int fillExtent[6];    // extents of fill
  myTree->fill->GetExtent(fillExtent);

  for(int i=0; i<3; i++){
    startPos[i] = spos[i];
    startDir[i] = sdir[i];
  }
  startRad = r;

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
  firstSeg = new veinSeg(this);
  lastSeg = firstSeg;

  // update length
  curLength += firstSeg->length;

  if (firstSeg->length == 0.0){
    failSeg = true;
  }

  int breastExtent[6];
  myTree->breast->GetExtent(breastExtent);

  // check if at ROI boundary by seeing if any neighboring voxels are outside ROI
  double* thePos = lastSeg->endPos;
  double pcoords[3];
  int inVol;
  inVol = myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
  if(inVol){
    // check if near phantom boundary
    if(invox[0] <= breastExtent[0] || invox[0] >= breastExtent[1] ||
       invox[1] <= breastExtent[2] || invox[1] >= breastExtent[3] ||
       invox[2] <= breastExtent[4] || invox[2] >= breastExtent[5]){
      failSeg = true;
    }
  } else {
    failSeg = true;
  }
  if(!failSeg){
    for(int a=-1; a<=1; a++){
      for(int b=-1; b<=1; b++){
	for(int c=-1; c<=1; c++){
	  unsigned char* p =
	    static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	  if(p[0] == myTree->tissue->skin || p[0] == myTree->tissue->bg){
	    edgeSeg = true;
	  }
	}
      }
    }
  }

  // generate more segments until proper length
  while(curLength < length && !failSeg && !edgeSeg){
    lastSeg->nextSeg = new veinSeg(lastSeg);  // the lastSeg in parenthesis is used to fill variables including prevSeg ptr
    lastSeg = lastSeg->nextSeg;
    curLength += lastSeg->length;
    if(lastSeg->length == 0.0){
      failSeg = true;
    }
    // check if at ROI boundary
    thePos = lastSeg->endPos;
    inVol = myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
    if(inVol){
      if(invox[0] <= breastExtent[0] || invox[0] >=breastExtent[1]||
	 invox[1] <= breastExtent[2] || invox[1] >=breastExtent[3] ||
	 invox[2] <= breastExtent[4] || invox[2] >=breastExtent[5]){
	failSeg = true;
      }
    } else {
      failSeg = true;
    }
    if(!failSeg){
      for(int a=-1; a<=1; a++){
	for(int b=-1; b<=1; b++){
	  for(int c=-1; c<=1; c++){
	    unsigned char* p =
	      static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	    if(p[0] == myTree->tissue->skin || p[0] == myTree->tissue->bg){
	      edgeSeg = true;
	    }
	  }
	}
      }
    }
  }

  // insert segments into phantom and update fill map
  veinSeg* mySeg = firstSeg;
  veinSeg* prevSeg;
  do{
    mySeg->updateMap();
#pragma omp parallel for collapse(3)
    for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
      for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
        for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
          double* v = static_cast<double*>(myTree->fill->GetScalarPointer(a,b,c));
          if(v[0] > 0.0){
            double dist;
            // voxel in ROI
            // voxel location
            vtkIdType id;
            int coord[3];
            coord[0] = a;
            coord[1] = b;
            coord[2] = c;
            id = myTree->fill->ComputePointId(coord);
            // get spatial coordinates of point
            double pos[3];
            myTree->fill->GetPoint(id,pos);
            dist = vtkMath::Distance2BetweenPoints(mySeg->endPos, pos);
            if(dist < v[0]){
              // update minimum distance
              v[0] = dist;
            }
          }
        }
      }
    }
    prevSeg = mySeg;
    mySeg = mySeg->nextSeg;
  } while(prevSeg != mySeg);

  // fill in end of branch variables
  for(int i=0; i<3; i++){
    endPos[i] = lastSeg->endPos[i];
    endDir[i] = lastSeg->endDir[i];
  }
  endRad = lastSeg->endRad;
  
  length = curLength;

  // set number of children and generate them
  nChild = setChild();

  // debug
  //cout << "br " << id << " len = " << length << " lev = " << level << " gen = " << gen <<  " child = " << nChild << "\n";
  //cout << "start dir = " << startDir[0] << "," << startDir[1] << "," << startDir[2] << "\tstart pos = " << startPos[0]<< "," << startPos[1] << "," <<startPos[2] << "\n";
  //cout << "end dir = " << endDir[0]<< "," << endDir[1] << "," <<endDir[2] << "\tend pos = " << endPos[0]<< "," << endPos[1] << "," <<endPos[2] << "\n";

  if(failSeg){
    nChild = 0;
    //std::cout << "Segment generation failure for branch" << id << std::endl;
  }

  if(edgeSeg){
    nChild = 0;
    std::cout << "ROI edge collision for branch " << id << std::endl;
  }

  if (nChild == 0){
    firstChild = nullptr;
    secondChild = nullptr;
  } else {
    // bifurcate
    // pick radii
    double radii[2];
    double thetas[2];
    setRadiiThetas(radii,thetas);
    // setup first child with level equal to current level
    firstChild = new veinBr(this,level,gen+1,radii[0],thetas[0]);
    firstChild->sibBranch = new veinBr(this,firstChild,level+1,gen+1,radii[1],thetas[1]);
    secondChild = firstChild->sibBranch;
  }
}

// constructor for first child branch of a parent branch
veinBr::veinBr(veinBr* par, unsigned int lev, unsigned int g, double r, double theta){

  int maxSegTry = 10;
  int numSegTry = 0;

  bool failSeg = false;
  bool edgeSeg = false;
  int invox[3];

  parent = par;
  sibBranch = nullptr;
  myTree = parent->myTree;

  int fillExtent[6];    // extents of fill
  myTree->fill->GetExtent(fillExtent);

  for(int i=0; i<3; i++){
    startPos[i] = parent->endPos[i];
  }
  startRad = r;
  level = lev;
  gen = g;
  id = myTree->numBranch;
  myTree->numBranch += 1;

  length = setLength();
  bool segSuccess = false;

  int breastExtent[6];
  myTree->breast->GetExtent(breastExtent);

  do{
    numSegTry++;
    setDir(startDir, theta);
    curLength = 0.0;

    // generate segments to fill branch
    firstSeg = new veinSeg(this);
    lastSeg = firstSeg;
  
    // update length
    curLength += firstSeg->length;

    if (firstSeg->length == 0.0){
      failSeg = true;
    }

    // check if at ROI boundary
    double* thePos = lastSeg->endPos;
    double pcoords[3];
    int inVol;
    inVol = myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
    if(inVol){
      if(invox[0] <= breastExtent[0] || invox[0] >= breastExtent[1] ||
	 invox[1] <= breastExtent[2] || invox[1] >= breastExtent[3] ||
	 invox[2] <= breastExtent[4] || invox[2] >= breastExtent[5]){
	failSeg = true;
      }
    } else {
      failSeg = true;
    }
    if(!failSeg){
      for(int a=-1; a<=1; a++){
	for(int b=-1; b<=1; b++){
	  for(int c=-1; c<=1; c++){
	    unsigned char* p =
	      static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	    if(p[0] == myTree->tissue->skin || p[0] == myTree->tissue->bg || p[0] == myTree->tissue->muscle){
	      edgeSeg = true;
	    }
	  }
	}
      }
    }

    // generate more segments until proper length
    while(curLength < length && !failSeg && !edgeSeg){
      lastSeg->nextSeg = new veinSeg(lastSeg);  // the lastSeg in parenthesis is used to fill variables including prevSeg ptr
      lastSeg = lastSeg->nextSeg;
      curLength += lastSeg->length;
      if(lastSeg->length == 0.0){
	failSeg = true;
      }
      // check if at ROI boundary
      thePos = lastSeg->endPos;
      inVol = myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
      if(inVol){
	if(invox[0] <= breastExtent[0] || invox[0] >= breastExtent[1] ||
	   invox[1] <= breastExtent[2] || invox[1] >= breastExtent[3] ||
	   invox[2] <= breastExtent[4] || invox[2] >= breastExtent[5]){
	  failSeg = true;
	}
      } else {
	failSeg = true;
      }
      if(!failSeg){
	for(int a=-1; a<=1; a++){
	  for(int b=-1; b<=1; b++){
	    for(int c=-1; c<=1; c++){
	      unsigned char* p =
		static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	      if(p[0] == myTree->tissue->skin || p[0] == myTree->tissue->bg || p[0] == myTree->tissue->muscle){
		edgeSeg = true;
	      }
	    }
	  }
	}
      }
    }

    if(failSeg || edgeSeg){
      //cout << "failed\n";
      // delete current segments and try again
      veinSeg* delSeg;
      while(firstSeg != lastSeg){
	delSeg = firstSeg;
	firstSeg = firstSeg->nextSeg;
	delete(delSeg);
      }
      delete(firstSeg);
    } else {
      //cout << "good!\n";
      segSuccess = true;
    }
  } while(!segSuccess && numSegTry < maxSegTry);

  // insert segments into phantom and update fill map
  veinSeg* mySeg = firstSeg;
  veinSeg* prevSeg;

  do{
    mySeg->updateMap();
#pragma omp parallel for collapse(3)
    for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
      for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
        for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
          double* v = static_cast<double*>(myTree->fill->GetScalarPointer(a,b,c));
          if(v[0] > 0.0){
            double dist;
            // voxel in ROI
            // voxel location
            vtkIdType id;
            int coord[3];
            coord[0] = a;
            coord[1] = b;
            coord[2] = c;
            id = myTree->fill->ComputePointId(coord);
            // get spatial coordinates of point
            double pos[3];
            myTree->fill->GetPoint(id,pos);
            dist = vtkMath::Distance2BetweenPoints(mySeg->endPos, pos);
            if(dist < v[0]){
              // update minimum distance
              v[0] = dist;
            }
          }
        }
      }
    }

    prevSeg = mySeg;
    mySeg = mySeg->nextSeg;
  } while(prevSeg != mySeg);
      
  // fill in end of branch variables
  for(int i=0; i<3; i++){
    endPos[i] = lastSeg->endPos[i];
    endDir[i] = lastSeg->endDir[i];
  }
  endRad = lastSeg->endRad;

  length = curLength;

  // set number of children and generate them
  nChild = setChild();

  // debug
  //cout << "br " << id << " len = " << length << " lev = " << level << " gen = " << gen <<  " child = " << nChild << "\n";
  //cout << "start dir = " << startDir[0]<< "," << startDir[1] << "," <<startDir[2] << "\tstart pos = " << startPos[0]<< "," << startPos[1] << "," <<startPos[2] << "\n";
  //cout << "end dir = " << endDir[0]<< "," << endDir[1] << "," <<endDir[2] << "\tend pos = " << endPos[0]<< "," << endPos[1] << "," <<endPos[2] << "\n";

  //if(failSeg){
  //  nChild = 0;
  //  std::cout << "Segment generation failure for branch" << id << std::endl;
  //}

  //if(edgeSeg){
  //  nChild = 0;
  //  std::cout << "ROI edge collision for branch " << id << std::endl;
  //}

  if(!segSuccess){
    nChild = 0;
    //std::cout << "Segment generation failure for branch " << id<< std::endl;
  }
  
  if (nChild == 0){
    firstChild = nullptr;
    secondChild = nullptr;
  } else {
    // pick radii
    double radii[2];
    double thetas[2];
    setRadiiThetas(radii,thetas);
		
    // setup first child with level equal to current level
    firstChild = new veinBr(this,level,gen+1,radii[0],thetas[0]);
    firstChild->sibBranch = new veinBr(this,level+1,gen+1,radii[1],radii[1]);
    secondChild = firstChild->sibBranch;

  }
}

// constructor for subsequent children (not first child) of a parent branch
veinBr::veinBr(veinBr* par, veinBr* par2, unsigned int lev, unsigned int g, double r, double theta){

  int maxSegTry = 10;
  int numSegTry = 0;

  bool failSeg = false;
  bool edgeSeg = false;

  int invox[3];

  parent = par;
  sibBranch = par2;

  myTree = parent->myTree;

  int fillExtent[6];    // extents of fill
  myTree->fill->GetExtent(fillExtent);

  for(int i=0; i<3; i++){
    startPos[i] = parent->endPos[i];
  }
  startRad = r;
  level = lev;
  gen = g;
  id = myTree->numBranch;
  myTree->numBranch += 1;
  
  length = setLength();
  bool segSuccess = false;

  int breastExtent[6];
  myTree->breast->GetExtent(breastExtent);

  do{
    numSegTry++;
    setDir(startDir, theta);
    curLength = 0.0;

    // generate segments to fill branch
    firstSeg = new veinSeg(this);
    lastSeg = firstSeg;

    // update length
    curLength += firstSeg->length;
  
    if (firstSeg->length == 0.0){
      failSeg = true;
    }

    // check if at ROI boundary
    double* thePos = lastSeg->endPos;
    double pcoords[3];
    int inVol;
    inVol = myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
    if(inVol){
      if(invox[0] <= breastExtent[0] || invox[0] >= breastExtent[1] ||
	 invox[1] <= breastExtent[2] || invox[1] >= breastExtent[3] ||
	 invox[2] <= breastExtent[4] || invox[2] >= breastExtent[5]){
	failSeg = true;
      }
    } else {
      failSeg = true;
    }
    if(!failSeg){
      for(int a=-1; a<=1; a++){
	for(int b=-1; b<=1; b++){
	  for(int c=-1; c<=1; c++){
	    unsigned char* p =
	      static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	    if(p[0] == myTree->tissue->skin || p[0] == myTree->tissue->bg || p[0] == myTree->tissue->muscle){
	      edgeSeg = true;
	    }
	  }
	}
      }
    }

    // generate more segments until proper length
    while(curLength < length && !failSeg && !edgeSeg){
      lastSeg->nextSeg = new veinSeg(lastSeg);  // the lastSeg in parenthesis is used to fill variables including prevSeg ptr
      lastSeg = lastSeg->nextSeg;
      curLength += lastSeg->length;
      if(lastSeg->length == 0.0){
	failSeg = true;
      }
      // check if at ROI boundary
      thePos = lastSeg->endPos;
      inVol = myTree->breast->ComputeStructuredCoordinates(thePos, invox, pcoords);
      if(inVol){
	if(invox[0] <= breastExtent[0] || invox[0] >= breastExtent[1] ||
	   invox[1] <= breastExtent[2] || invox[1] >= breastExtent[3] ||
	   invox[2] <= breastExtent[4] || invox[2] >= breastExtent[5]){
	  failSeg = true;
	}
      } else {
	failSeg = true;
      }
      if(!failSeg){
	for(int a=-1; a<=1; a++){
	  for(int b=-1; b<=1; b++){
	    for(int c=-1; c<=1; c++){
	      unsigned char* p =
		static_cast<unsigned char*>(myTree->breast->GetScalarPointer(invox[0]+a,invox[1]+b,invox[2]+c));
	      if(p[0] == myTree->tissue->skin || p[0] == myTree->tissue->bg || p[0] == myTree->tissue->muscle){
		edgeSeg = true;
	      }
	    }
	  }
	}
      }
    }

    if(failSeg || edgeSeg){
      //cout << "failed\n";
      // delete current segments and try again
      veinSeg* delSeg;
      while(firstSeg != lastSeg){
        delSeg = firstSeg;
        firstSeg = firstSeg->nextSeg;
        delete(delSeg);
      }
      delete(firstSeg);
    } else {
      //cout << "good!\n";
      segSuccess = true;
    }
  } while(!segSuccess && numSegTry < maxSegTry);

  // insert segments into phantom and update fill map
  veinSeg* mySeg = firstSeg;
  veinSeg* prevSeg;

  do{
    mySeg->updateMap();
    // update density map
#pragma omp parallel for collapse(3)
    for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
      for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
        for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
          double* v = static_cast<double*>(myTree->fill->GetScalarPointer(a,b,c));
          if(v[0] > 0.0){
            double dist;
            // voxel in ROI
            // voxel location
            vtkIdType id;
            int coord[3];
            coord[0] = a;
            coord[1] = b;
            coord[2] = c;
            id = myTree->fill->ComputePointId(coord);
            // get spatial coordinates of point
            double pos[3];
            myTree->fill->GetPoint(id,pos);
            dist = vtkMath::Distance2BetweenPoints(mySeg->endPos, pos);
            if(dist < v[0]){
              // update minimum distance
	      v[0] = dist;
            }
          }
        }
      }
    }

    prevSeg = mySeg;
    mySeg = mySeg->nextSeg;
  } while(prevSeg != mySeg);

  // fill in end of branch variables
  for(int i=0; i<3; i++){
    endPos[i] = lastSeg->endPos[i];
    endDir[i] = lastSeg->endDir[i];
  }
  endRad = lastSeg->endRad;
  
  length = curLength;

  // set number of children and generate them
  nChild = setChild();
  
  // debug
  //cout << "br " << id << " len = " << length << " lev = " << level << " gen = " << gen <<  " child = " << nChild << "\n";
  //cout << "start dir = " << startDir[0]<< "," << startDir[1] << "," <<startDir[2] << "\tstart pos = " << startPos[0]<< "," << startPos[1] << "," <<startPos[2] << "\n";
  //cout << "end dir = " << endDir[0]<< "," << endDir[1] << "," <<endDir[2] << "\tend pos = " << endPos[0]<< "," << endPos[1] << "," <<endPos[2] << "\n";

  if(!segSuccess){
    nChild = 0;
    //std::cout << "Segment generation failure for branch " << id << std::endl;
  }

  if (nChild == 0){
    firstChild = nullptr;
    secondChild = nullptr;
  } else {
    // pick radii and angles
    double radii[2];
    double thetas[2];
    setRadiiThetas(radii, thetas);

    // setup first child with level equal to current level
    firstChild = new veinBr(this,level,gen+1,radii[0],thetas[0]);
    firstChild->sibBranch = new veinBr(this,level+1,gen+1,radii[1],thetas[1]);
    secondChild = firstChild->sibBranch;

  }
}

double veinBr::setLength(void){
  double len;
  double randVal = myTree->u01();
  double baseLen = myTree->baseLength;
  double lenShrink = myTree->opt["vesselBr.lenShrink"].as<double>();
  double lenRange = myTree->opt["vesselBr.lenRange"].as<double>();

  len = baseLen*pow(lenShrink,level);

  len = len - lenRange*len + randVal*2*lenRange*len;
  return(len);
}

unsigned int veinBr::setChild(void){
  // determine number of child branches
  // if small enough, no children
  double minRad = myTree->opt["vesselBr.childMinRad"].as<double>();
  if(endRad < minRad){
    return(0);
  }

  // check for max number branches
  if(myTree->numBranch >= myTree->maxBranch){
    return(0);
  }

  // define maximum generation
  unsigned int maxGen = myTree->opt["vesselTree.maxGen"].as<uint>();
  if(gen > maxGen){
    return(0);
  }

  // default to 2 children
  return(2);
}

void veinBr::setRadiiThetas(double* radii, double* thetas){
  // set radii of child branches based on radii of the parent

  double minFrac = myTree->opt["vesselBr.minRadFrac"].as<double>();
  double maxFrac = myTree->opt["vesselBr.maxRadFrac"].as<double>();

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

void veinBr::setDir(double* sdir, double theta){

  const double pi = boost::math::constants::pi<double>();

  // set initial direction of child branches
  double dir[3];
  double tempV[3];
  double basis1[3];
  double basis2[3];

  double rotateJitter = myTree->opt["ductBr.rotateJitter"].as<double>();
  double rotate;
  double minAngleSep = 0.1;

  if(sibBranch == nullptr){
    // this is the first child
    // random rotation about parent direction
    rotate = 2*pi*myTree->u01();
    azimuth = rotate;
  } else {
    // this is the second child, keep away from sibling
    do{
      rotate = 2*pi*myTree->u01();
    } while(fabs(rotate-sibBranch->azimuth) < minAngleSep);
    //rotate = rotate - rotateJitter + randVal*2*rotateJitter;
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
veinBr::~veinBr(){
  // delete segments
  veinSeg* delSeg;
  while(firstSeg != lastSeg){
    delSeg = firstSeg;
    firstSeg = firstSeg->nextSeg;
    delete(delSeg);
  }
  delete(firstSeg);

  // delete child branches
  veinBr* delBranch;
  while(firstChild != secondChild){
    delBranch = firstChild;
    firstChild = firstChild->secondChild;
    delete(delBranch);
  }
  delete(firstChild);
}

// constructor for first segment
veinSeg::veinSeg(veinBr* br){
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
veinSeg::veinSeg(veinSeg* pr){
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

void veinSeg::makeSeg(){

  const double pi = boost::math::constants::pi<double>();

  double segFrac = myBranch->myTree->opt["ductSeg.segFrac"].as<double>();
  unsigned int numTry = myBranch->myTree->opt["vesselSeg.numTry"].as<uint>();
  unsigned int maxTry = myBranch->myTree->opt["vesselSeg.maxTry"].as<uint>();
  unsigned int absMaxTry = myBranch->myTree->opt["vesselSeg.absMaxTry"].as<uint>();
  double maxRad = myBranch->myTree->opt["vesselSeg.maxCurvRad"].as<double>();
  // debug
  maxRad = maxRad/(myBranch->level+1.0);
  double angleMax =  pi*myBranch->myTree->opt["vesselSeg.maxCurvFrac"].as<double>();
  double roiStep = myBranch->myTree->opt["vesselSeg.roiStep"].as<double>();
  double densityWt = myBranch->myTree->opt["vesselSeg.densityWt"].as<double>();
  double angleWt = myBranch->myTree->opt["vesselSeg.angleWt"].as<double>();
  double dirWt = myBranch->myTree->opt["vesselSeg.dirWt"].as<double>();
  double prefDir[3]; // preferential direction of growth
  for(int i=0; i<3; i++){
    prefDir[i] = myBranch->myTree->nipplePos[i] - startPos[i];
  }
  vtkMath::Normalize(prefDir);
  double maxEndRad = myBranch->myTree->opt["vesselSeg.maxEndRad"].as<double>();
  double minEndRad = myBranch->myTree->opt["vesselSeg.minEndRad"].as<double>();

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
  int inVol;        // Structured coordinates check

  double travelDist;
  double travelStep = 1.0;

  // breast FOV for checking if segment in FOV
  double breastOrigin[3];
  double breastSpacing[3];
  int breastDim[3];
  double breastFOV[6];
  int breastExtent[6];

  myBranch->myTree->breast->GetOrigin(breastOrigin);
  myBranch->myTree->breast->GetSpacing(breastSpacing);
  //myBranch->myTree->breast->GetDimensions(breastDim);
  myBranch->myTree->breast->GetExtent(breastExtent);

  for(int i=0; i<3; i++){
    breastFOV[2*i] = breastOrigin[i]+(double)breastExtent[2*i]*breastSpacing[i]-breastSpacing[i]/2.0;
    breastFOV[2*i+1] = breastOrigin[i]+(double)breastExtent[2*i+1]*breastSpacing[i]+breastSpacing[i]/2.0;
  }

  double bestCurv[3];	// best center of curvature found so far
  double bestRadius;
  double bestCost;  // best cost found so far
  double cost;		// current cost

  double density;	// change in density of vein structure due to new segment

  // determine proposed segment length
  length = myBranch->length*segFrac;

  if(length > (myBranch->length - myBranch->curLength)){
    // truncate
    length = myBranch->length - myBranch->curLength;
  }

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
	    inVol = myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, myVoxel, pcoords);
	    if(inVol){
	      unsigned char* p =
		static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(myVoxel[0],myVoxel[1],myVoxel[2]));
	      bool inBreast = true;
	      if(p[0] == myBranch->myTree->tissue->skin || p[0] == myBranch->myTree->tissue->bg){
		inBreast = false;
	      }
	      if(!inBreast){
		inROI = false;
	      }
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
	  inVol = myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, myVoxel, pcoords);
	  if(inVol){
	    unsigned char* p =
	      static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(myVoxel[0],myVoxel[1],myVoxel[2]));
	    bool inBreast = true;
	    if(p[0] == myBranch->myTree->tissue->skin || p[0] == myBranch->myTree->tissue->bg){
	      inBreast = false;
	    }
	    if(!inBreast){
	      inROI = false;
	    }
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

	  // reduction in squared distance to arteries in ROI
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
	  
	  // test if heading for edge
          travelDist = 0.0;
	  bool inBreast = true;

          while(inBreast){
            double currPos[3];
            for(int i=0; i<3; i++){
              currPos[i] = checkPos[i] + travelDist*endDir[i];
            }
            inVol = myBranch->myTree->breast->ComputeStructuredCoordinates(currPos, myVoxel, pcoords);
            if(inVol){
              unsigned char* p =
                static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(myVoxel[0],myVoxel[1],myVoxel[2]));
	      if(p[0] == myBranch->myTree->tissue->skin || p[0] == myBranch->myTree->tissue->bg){
		inBreast = false;
	      } else {
		travelDist += travelStep;
	      } 
            } else {
	      inBreast = false;
            }
          }

	  // prefDir towards nipple
	  for(int i=0; i<3; i++){
	    prefDir[i] = myBranch->myTree->nipplePos[i] - checkPos[i];
	  }
	  vtkMath::Normalize(prefDir);

	  bestCost = densityWt*density - angleWt*vtkMath::Dot(endDir,prefDir) - dirWt*travelDist;
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

	  // test if heading for edge
	  double travelDist = 0.0;
	  // step size
	  double travelStep = 1.0;
	  bool inBreast = true;
	  
	  while(inBreast){
	    double currPos[3];
	    for(int i=0; i<3; i++){
	      currPos[i] = checkPos[i] + travelDist*endDir[i];
	    }
	    inVol = myBranch->myTree->breast->ComputeStructuredCoordinates(currPos, myVoxel, pcoords);
	    if(inVol){
	      unsigned char* p =
		static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(myVoxel[0],myVoxel[1],myVoxel[2]));
	      if(p[0] == myBranch->myTree->tissue->skin || p[0] == myBranch->myTree->tissue->bg){
		inBreast = false;
	      } else {
		travelDist += travelStep;
	      }
	    } else {
	      inBreast = false;
	    }
	  }

	  // prefDir towards nipple
	  for(int i=0; i<3; i++){
            prefDir[i] = myBranch->myTree->nipplePos[i] - checkPos[i];
          }
	  vtkMath::Normalize(prefDir);

	  cost = densityWt*density - angleWt*vtkMath::Dot(endDir,prefDir) - dirWt*travelDist;
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
      // debug
      //cout << "Reducing segment length\n";
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
    // debug
    //cout << "Cost= " << bestCost << "\n";
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
    //updateMap();
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

void veinSeg::setShape(){
  // cubic polynomial v(0)d^3+v(1)d^2+v(2)d+v(3)
  // d = [0,length)
  shape[3] = startRad;
  shape[2] = startDeriv;
  double c = endRad-shape[2]*length-shape[3];
  shape[0] = (endDeriv-shape[2]-2*c/length)/length/length;
  shape[1] = c/length/length-shape[0]*length;
}

double veinSeg::getRadius(double t){
  // return segment radius at position t (mm)
  if(t>=0 && t<=length){
    return(shape[0]*t*t*t+shape[1]*t*t+shape[2]*t+shape[3]);
  } else {
    return(0);
  }
}

void veinSeg::updateMap(){

  const double pi = boost::math::constants::pi<double>();

  // update voxelized map of arteries
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

  // position along segment, radius and angle
  //double lpos=0.0;

  // calculate number of for loops for openMP
  int lIter = (int)(ceil(length/ls));
	
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
      int inVol;
			
      if(rpos < step){
	// only check current voxel assume angle = 0
	for(int i=0; i<3; i++){
	  checkPos[i] = currentPos[i] + rpos*(-1*cos(0.0)*lbasis2[i] + sin(0.0)*basis3[i]);
	}
	inVol = myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, checkIdx, pcoords);
	if(inVol){
	  // set voxel to vein
	  unsigned char* p =
	    static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(checkIdx[0],checkIdx[1],checkIdx[2]));
	  p[0] = myBranch->myTree->tissue->vein;
	}
      } else {
	// angle step
	double as = step/rpos;
	double apos = 0.0;
	while(apos < 2*pi){
	  for(int i=0; i<3; i++){
	    checkPos[i] = currentPos[i] + rpos*(-1*cos(apos)*lbasis2[i] + sin(apos)*basis3[i]);
	  }
	  inVol = myBranch->myTree->breast->ComputeStructuredCoordinates(checkPos, checkIdx, pcoords);
	  if(inVol){
	    // set voxel to vein
	    unsigned char* p =
	      static_cast<unsigned char*>(myBranch->myTree->breast->GetScalarPointer(checkIdx[0],checkIdx[1],checkIdx[2]));
	    p[0] = myBranch->myTree->tissue->vein;
	  }
	  apos += as;
	}
      }
      rpos += rs;
    }  
    //lpos += ls;
  }
}
