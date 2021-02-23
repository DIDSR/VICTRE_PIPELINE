/*
 * create_artery.cxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

// function to generate an arterial tree

#include "createArtery.hxx"

using namespace std;
namespace po = boost::program_options;

/* This function creates arterial network, inserts it into the segmented
 * breast and saves the tree */
void generate_artery(vtkImageData* breast, po::variables_map vm, int* boundBox,
		     tissueStruct* tissue, double* sposPtr, double* sdirPtr, double* nipplePos, int seed, int mainSeed, bool firstTree){

  char arteryFilename[128];
  std::string outputDir = vm["base.outputDir"].as<std::string>();
  sprintf(arteryFilename,"%s/p_%d_arteryFill.vti", outputDir.c_str(), mainSeed);
  
  double spos[3];
  double sdir[3];

  for(int i=0; i<3; i++){
    spos[i] = sposPtr[i];
    sdir[i] = sdirPtr[i];
  }

  // declare arteryTreeInit struct and fill information
  arteryTreeInit treeInit;

  treeInit.seed = seed;

  // bounds of artery simulation derived from breast structure
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

  treeInit.nFill[0] = vm["vesselTree.nFillX"].as<uint>();
  treeInit.nFill[1] = vm["vesselTree.nFillY"].as<uint>();
  treeInit.nFill[2] = vm["vesselTree.nFillZ"].as<uint>();

  for(int i=0; i<3; i++){
    treeInit.nipplePos[i] = nipplePos[i];
  }

  treeInit.boundBox = boundBox;

  treeInit.tissue = tissue;

  treeInit.breast = breast;

  // create arterial tree
  arteryTree myTree(vm, &treeInit);

  // root of tree
  double srad = vm["vesselTree.initRad"].as<double>();

  // initialize fill map based on distance to start position if first tree, else load current fill

  if(firstTree){
    // initialize fill map based on distance to start position
    int fillExtent[6];
    myTree.fill->GetExtent(fillExtent);
    for(int a=fillExtent[0]; a<=fillExtent[1]; a++){
      for(int b=fillExtent[2]; b<=fillExtent[3]; b++){
	for(int c=fillExtent[4]; c<=fillExtent[5]; c++){
	  double* v = static_cast<double*>(myTree.fill->GetScalarPointer(a,b,c));
	  // set distance to 0 if fill voxel is not in breast, otherwise
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
	  bool inBreast = true;
	  unsigned char voxelVal = breastVal[breast->FindPoint(pos)];
	  if(voxelVal == tissue->skin || voxelVal == tissue->bg){
	    inBreast = false;
	  }
	  if(inBreast){
	    // inside breast
	    v[0] = vtkMath::Distance2BetweenPoints(spos,pos);
	  } else {
	    // outside breast, set distance to zero
	    v[0] = 0.0;
	  }
	}
      }
    }
  } else {
    // debug
    vtkSmartPointer<vtkXMLImageDataReader> fillReader =
      vtkSmartPointer<vtkXMLImageDataReader>::New();
    fillReader->SetFileName(arteryFilename);
    fillReader->Update();
    vtkSmartPointer<vtkImageData> holder = 
      vtkSmartPointer<vtkImageData>::New();
    holder = fillReader->GetOutput();
    myTree.fill->DeepCopy(holder);
		
    // load from file
    //vtkSmartPointer<vtkXMLImageDataReader> fillReader =
    //		vtkSmartPointer<vtkXMLImageDataReader>::New();
    //fillReader->SetFileName("arteryFill.vti");
    //fillReader->Update();
    //myTree.fill = fillReader->GetOutput();
		
  }

  myTree.head = new arteryBr(spos, sdir, srad, &myTree);

  // save density map
  vtkSmartPointer<vtkXMLImageDataWriter> fillWriter =
    vtkSmartPointer<vtkXMLImageDataWriter>::New();

  fillWriter->SetFileName(arteryFilename);
#if VTK_MAJOR_VERSION <= 5
  fillWriter->SetInput(myTree.fill);
#else
  fillWriter->SetInputData(myTree.fill);
#endif
  fillWriter->Write();

  //  cout << "Finished artery " << myTree.id << ", " <<
  //  myTree.numBranch << " branches created\n";

  return;
}















