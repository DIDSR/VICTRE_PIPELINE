/*
 * create_duct.cxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#include "createDuct.hxx"

using namespace std;
namespace po = boost::program_options;

/* This function creates a duct tree within a given compartment, inserts it into the segmented
 * breast and saves the tree */

void generate_duct(vtkImageData* breast, po::variables_map vm, vtkPoints* TDLUloc, vtkDoubleArray* TDLUattr, 
		   unsigned char compartmentId, int* boundBox, tissueStruct* tissue, double* sposPtr, double* sdirPtr, int seed){

  double spos[3];
  double sdir[3];

  for(int i=0; i<3; i++){
    spos[i] = sposPtr[i];
    sdir[i] = sdirPtr[i];
  }

  // declare ductTreeInit struct and fill information
  ductTreeInit treeInit;

  treeInit.seed = seed;

  // bounds of duct simulation derived from breast structure
  int startInd[3] = {boundBox[0], boundBox[2], boundBox[4]};
  int endInd[3] = {boundBox[1], boundBox[3], boundBox[5]};

  //startPos
  breast->GetPoint(breast->ComputePointId(startInd), treeInit.startPos);
  //endPos
  breast->GetPoint(breast->ComputePointId(endInd), treeInit.endPos);

  // size of voxels
  treeInit.nVox[0] = boundBox[1]-boundBox[0];
  treeInit.nVox[1] = boundBox[3]-boundBox[2];
  treeInit.nVox[2] = boundBox[5]-boundBox[4];

  treeInit.nFill[0] = vm["ductTree.nFillX"].as<uint>();
  treeInit.nFill[1] = vm["ductTree.nFillY"].as<uint>();
  treeInit.nFill[2] = vm["ductTree.nFillZ"].as<uint>();

  for(int i=0; i<3; i++){
    treeInit.prefDir[i] = sdir[i];
  }

  treeInit.boundBox = boundBox;

  treeInit.compartmentId = compartmentId;
  
  treeInit.tissue = tissue;
  
  treeInit.breast = breast;
	
  treeInit.TDLUloc = TDLUloc;
	
  treeInit.TDLUattr = TDLUattr;

  ductTree myTree(vm, &treeInit);

  // root of tree
  double srad = vm["ductTree.initRad"].as<double>();

  // initialize fill map based on distance to start position
  int fillExtent[6];
  myTree.fill->GetExtent(fillExtent);
#pragma omp parallel for
  for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
    for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
      for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
	double* v = static_cast<double*>(myTree.fill->GetScalarPointer(a,b,c));
	// set distance to 0 if fill voxel is not in compartment, otherwise
	// initialize with squared distance to tree base

	// fill voxel location
	vtkIdType id;
	int coord[3];
	coord[0] = a;
	coord[1] = b;
	coord[2] = c;
	id = myTree.fill->ComputePointId(coord);
	// get spatial coordinates of fill voxel
	double pos[3];
	myTree.fill->GetPoint(id,pos);
	// compare to nearest breast voxel id
	unsigned char* breastVal = static_cast<unsigned char *>(breast->GetScalarPointer());
	if(breastVal[breast->FindPoint(pos)] == compartmentId){
	  // inside compartment
	  v[0] = vtkMath::Distance2BetweenPoints(spos, pos);
	} else {
	  // outside compartment, set distance to zero
	  v[0] = 0.0;
	}
      }
    }
  }

  myTree.head = new ductBr(spos, sdir, srad, &myTree);

  //cout << "Finished duct tree " << myTree.id << " with " << myTree.TDLUloc->GetNumberOfPoints() << " TDLUs, " <<
  //  myTree.numBranch << " branches created\n";
  
  return;
}















