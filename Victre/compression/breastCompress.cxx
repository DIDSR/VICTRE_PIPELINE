/*! \file breastCompress.cxx
 *  \brief breastCompress main file
 *  \author Christian G. Graff
 *  \version 1.0
 *  \date 2018
 *  
 *  \copyright To the extent possible under law, the author(s) have
 *  dedicated all copyright and related and neighboring rights to this
 *  software to the public domain worldwide. This software is
 *  distributed without any warranty.  You should have received a copy
 *  of the CC0 Public Domain Dedication along with this software.
 *  If not, see <http://creativecommons.org/publicdomain/zero/1.0/>.
 * 
 */

#include "breastCompress.hxx"

namespace po = boost::program_options;

int main(int argc, char *argv[])
{

  double pi = vtkMath::Pi();

  // breast tissue parameters
  // Young's modulus, Poisson ratio, density

  double fatModulus = 5.0; // kPa (equivalent to kg/(mm*s^2))
  double fatPoisson = 0.49;
  double fatDensity = 0.000001; // kg/mm^3

  double glandModulus = 15.0; // kPa
  double glandPoisson = 0.49;
  double glandDensity = 0.000001; // kg/mm^3

  double massModulus = 20.0; // kPa
  double massPoisson = 0.49;
  double massDensity = 0.00001; // kg/mm^3

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
  tissue.vein = 225;
  tissue.muscle = 40;
  tissue.mass = 200;
  tissue.calc = 250;
  tissue.paddle = 50;

  // command line parameters
  po::options_description opts("All options");
  opts.add_options()("seed,s", po::value<int>()->required(), "input breast seed")("febiopath,f", po::value<std::string>(), "FeBiO path")("angle,a", po::value<double>()->default_value(0.0), "compression paddle angle (degrees)")("thickness,t", po::value<double>()->required(), "final thickness (mm)")
      // 0 -> autoremesh, else number of meshes
      ("meshes,m", po::value<int>()->default_value(0), "Number of times to (re)mesh")("dir,d", po::value<std::string>()->default_value("."), "work directory")("tetvolume,v", po::value<double>()->default_value(25.0), "maximum tetrahedral volume")("decimateFrac,r", po::value<double>()->default_value(0.99), "surface mesh decimation fraction")("help,h", "print help message");

  po::variables_map vm;

  // get command line arguments
  po::store(parse_command_line(argc, argv, opts), vm);

  if (vm.count("help"))
  {
    cout << "\n";
    cout << "breastCompress - breast phantom compression simulation\n\n";
    cout << opts << "\n";
    return 1;
  }

  // determine FeBiO path
  char febioPath[256];

  if (vm.count("febiopath"))
  {
    // specified on command line
    std::string febio = vm["febiopath"].as<std::string>();
    strcpy(febioPath, febio.c_str());
  }
  else
  {
    // check environment variable
    const char *svar = getenv("FEBIO_PATH");
    if (svar == NULL)
    {
      cout << "Error: FeBiO path not specified\n";
      return (EXIT_FAILURE);
    }
    else
    {
      strcpy(febioPath, svar);
    }
  }

  if (!vm.count("seed") || !vm.count("thickness"))
  {
    cout << "Command line options missing\n";
    cout << opts << "\n";
    return (EXIT_FAILURE);
  }

  int seed = vm["seed"].as<int>();
  // compression angle converted to radians
  double angle = vm["angle"].as<double>();
  angle = angle * pi / 180.0;

  bool rotate = false;
  if (fabs(angle) > 1e-9)
  {
    rotate = true;
  }

  double thickness = vm["thickness"].as<double>();

  // check for nonsensical thickness
  if (thickness < 10.0 || thickness > 200.0)
  {
    cerr << "Compressed thickness should be between 10 and 200 mm.\n";
    return (1);
  }
  int numMeshes = vm["meshes"].as<int>();
  bool autoRemesh = (numMeshes == 0);

  // mesh parameters
  //double decimateFrac = 0.99;
  double decimateFrac = vm["decimateFrac"].as<double>();

  if (decimateFrac < 0.0 || decimateFrac >= 1.0)
  {
    cerr << "Invalid mesh decimation fraction.\n";
    return (1);
  }

  tetgenbehavior tetPar;
  tetPar.plc = 1;
  tetPar.quality = 1;
  tetPar.coarsen = 1;
  tetPar.fixedvolume = 1;
  //tetPar.maxvolume = 25.0;
  tetPar.maxvolume = vm["tetvolume"].as<double>();
  tetPar.optlevel = 10;
  //tetPar.nobisect = 1;
  tetPar.neighout = 2; // trigger calculation of boundary face parents
  //  tetPar.nomergevertex = 1;

  int paddleNh = 18;
  int paddleNl = 25;
  int paddleNw = 25;

  std::string workDir = vm["dir"].as<std::string>();

  bool done = false;

  double minimumFinalThicknessSoFar = 99999.9; // Dummy value.
  double prevTopPaddleCenterDerotated[3] = {0.0, 0.0, 0.0};
  double prevBotPaddleCenterDerotated[3] = {0.0, 0.0, 0.0};
  double origPaddleWidth = 0.0;
  double origPaddleLength = 0.0;
  double paddleDist = 0.0;

  int meshCount = 0;
  const int MAX_MESHES = 18;

  // main loop
  while ((autoRemesh && !done && (meshCount < MAX_MESHES)) || (!autoRemesh && (meshCount < numMeshes)))
  {

    char inputHdrFilename[256];
    char inputgzFilename[256];
    char inputLocFilename[256];

    vtkSmartPointer<vtkImageData> input =
        vtkSmartPointer<vtkImageData>::New();

    if (meshCount == 0)
    {
      sprintf(inputHdrFilename, "%s/p_%d.mhd", workDir.c_str(), seed);
      sprintf(inputLocFilename, "%s/p_%d.loc", workDir.c_str(), seed);
      sprintf(inputgzFilename, "%s/p_%d.raw.gz", workDir.c_str(), seed);
    }
    else
    {
      sprintf(inputHdrFilename, "%s/p_%d_%d.mhd", workDir.c_str(), seed, meshCount - 1);
      sprintf(inputLocFilename, "%s/p_%d_%d.loc", workDir.c_str(), seed, meshCount - 1);
      sprintf(inputgzFilename, "%s/p_%d_%d.raw.gz", workDir.c_str(), seed, meshCount - 1);
    }

    // parse header
    double inputSpacing[3];
    int inputDim[3];
    double inputOrigin[3];

    FILE *inputHdr = fopen(inputHdrFilename, "r");
    if (inputHdr == NULL)
    {
      cerr << "Unable to open header file " << inputHdrFilename << "\n";
      return (EXIT_FAILURE);
    }
    else
    {
      char hdrString[256];
      // dummy variable to receive fgets return val
      char *sReturn;
      for (int i = 0; i < 6; i++)
      {
        sReturn = fgets(hdrString, 256, inputHdr);
      }
      sReturn = fgets(hdrString, 256, inputHdr);
      sscanf(hdrString, "%*s %*s %lf %lf %lf", &inputOrigin[0], &inputOrigin[1], &inputOrigin[2]);
      sReturn = fgets(hdrString, 256, inputHdr);
      sReturn = fgets(hdrString, 256, inputHdr);
      sscanf(hdrString, "%*s %*s %lf %lf %lf", &inputSpacing[0], &inputSpacing[1], &inputSpacing[2]);
      sReturn = fgets(hdrString, 256, inputHdr);
      sscanf(hdrString, "%*s %*s %d %d %d", &inputDim[0], &inputDim[1], &inputDim[2]);
      fclose(inputHdr);
    }

    input->SetDimensions(inputDim);
    input->SetSpacing(inputSpacing);
    input->SetOrigin(inputOrigin);

    // allocate memory and load from gz file
    input->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

    gzFile gzf = gzopen(inputgzFilename, "rb");

    if (gzf == NULL)
    {
      cerr << "Unable to open gz data file" << inputgzFilename << "\n";
      return (EXIT_FAILURE);
    }
    else
    {
      unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer());
      unsigned char *buffer;
      int bytesRead;
      long long int totalBytesRead = 0;
      buffer = (unsigned char *)malloc(inputDim[1] * inputDim[2] * sizeof(char));
      gzbuffer(gzf, inputDim[1] * inputDim[2]);

      for (int i = 0; i < inputDim[0]; i++)
      {
        bytesRead = gzread(gzf, buffer, inputDim[1] * inputDim[2]);
        if (bytesRead < 0)
        {
          cout << gzerror(gzf, &bytesRead) << "\n";
        }
        else
        {
          totalBytesRead += bytesRead;
        }
        memcpy(p, buffer, inputDim[1] * inputDim[2]);
        p += inputDim[1] * inputDim[2];
      }

      gzclose(gzf);
      free(buffer);
      cout << "Total bytes read = " << totalBytesRead << "\n";
    }

    // create list of lesion insertion locations
    vtkSmartPointer<vtkPoints> lesionLoc =
        vtkSmartPointer<vtkPoints>::New();

    FILE *inLocFile = fopen(inputLocFilename, "rt");
    char line[128];
    while (fgets(line, 128, inLocFile) != NULL)
    {
      double pt[3];
      char *sp;
      sp = strtok(line, " ,");
      for (int i = 0; i < 3; i++)
      {
        pt[i] = atof(sp);
        sp = strtok(NULL, " ,");
      }
      lesionLoc->InsertNextPoint(pt);
    }
    fclose(inLocFile);

    int numLesionLoc = lesionLoc->GetNumberOfPoints();

    // create lower resolution mask for boundary meshing
    vtkSmartPointer<vtkImageData> mask =
        vtkSmartPointer<vtkImageData>::New();

    double voxelSize = 1.0;
    mask->SetSpacing(voxelSize, voxelSize, voxelSize);

    double inputBounds[6];
    input->GetBounds(inputBounds);

    //std::cout<<"\nBegining mask search of breast surface.\n\n";

    int maskDim[3];

    // pad mask around breast
    for (int i = 0; i < 3; i++)
    {
      maskDim[i] = (int)ceil((inputBounds[2 * i + 1] - inputBounds[2 * i]) / voxelSize) + (int)ceil(75.0 / voxelSize);
    }

    mask->SetExtent(0, maskDim[0] - 1, 0, maskDim[1] - 1, 0, maskDim[2] - 1);

    double maskOrigin[3];
    input->GetOrigin(maskOrigin);

    // Needed so that breast is completely enclosed on all sides by background
    for (int i = 0; i < 3; i++)
    {
      maskOrigin[i] -= 75.0 / 2;
    }

    mask->SetOrigin(maskOrigin);

    mask->AllocateScalars(VTK_UNSIGNED_CHAR, 1);
#pragma omp parallel for collapse(3)
    for (int k = 0; k < maskDim[2]; k++)
    {
      for (int j = 0; j < maskDim[1]; j++)
      {
        for (int i = 0; i < maskDim[0]; i++)
        {
          unsigned char *q = static_cast<unsigned char *>(mask->GetScalarPointer(i, j, k));
          bool foundIn = false;
          // check 27 points on voxel surface to check if any intersection with phantom
          for (int a = -1; a <= 1; a++)
          {
            for (int b = -1; b <= 1; b++)
            {
              for (int c = -1; c <= 1; c++)
              {
                double pt[3] = {maskOrigin[0] + i * voxelSize + a * voxelSize / 2,
                                maskOrigin[1] + j * voxelSize + b * voxelSize / 2,
                                maskOrigin[2] + k * voxelSize + c * voxelSize / 2};
                int abc[3];
                double pcoords[3];
                int inVol = input->ComputeStructuredCoordinates(pt, abc, pcoords);
                if (inVol)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(abc));
                  if (*p != tissue.bg && *p != tissue.paddle)
                  {
                    *q = 2;
                    foundIn = true;
                  }
                }
              }
            }
          }
          if (!foundIn)
          {
            *q = 0;
          }
        }
      }
    }

    // create bounding surface using discrete marching cubes
    vtkSmartPointer<vtkDiscreteMarchingCubes> surfaceMaker =
        vtkSmartPointer<vtkDiscreteMarchingCubes>::New();

    surfaceMaker->SetInputData(mask);
    //surfaceMaker->ComputeNormalsOn();
    surfaceMaker->SetNumberOfContours(1);
    surfaceMaker->SetValue(0, 2);
    surfaceMaker->Update();

    vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
        vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();

    smoother->SetInputConnection(surfaceMaker->GetOutputPort());
    smoother->SetNumberOfIterations(20);
    smoother->BoundarySmoothingOff();
    smoother->FeatureEdgeSmoothingOff();
    smoother->SetFeatureAngle(120.0);
    smoother->SetPassBand(0.001);
    smoother->NonManifoldSmoothingOn();
    smoother->NormalizeCoordinatesOn();
    smoother->Update();

    // clean mesh
    vtkSmartPointer<vtkCleanPolyData> cleanPoly =
        vtkSmartPointer<vtkCleanPolyData>::New();
    //cleanPoly->PointMergingOn();
    //cleanPoly->SetTolerance(0.002);
    cleanPoly->SetInputConnection(smoother->GetOutputPort());

    // decimate mesh
    vtkSmartPointer<vtkQuadricDecimation> decimate =
        vtkSmartPointer<vtkQuadricDecimation>::New();
    decimate->SetTargetReduction(decimateFrac);
    decimate->SetInputConnection(cleanPoly->GetOutputPort());
    decimate->Update();

    vtkSmartPointer<vtkPolyData> breastMeshPoly =
        vtkSmartPointer<vtkPolyData>::New();

    breastMeshPoly->ShallowCopy(decimate->GetOutput());

    /**********************
     * meshing with Tetgen
     *********************/

    std::cout << "\nMeshing with Tetgen\n\n";

    tetgenio tetin, tetout;

    tetin.firstnumber = 0;
    tetin.numberofpoints = breastMeshPoly->GetNumberOfPoints();

    tetin.pointlist = new REAL[tetin.numberofpoints * 3];

    // populate node list
    for (int i = 0; i < breastMeshPoly->GetNumberOfPoints(); i++)
    {
      double thisPoint[3];
      breastMeshPoly->GetPoint(i, thisPoint);
      tetin.pointlist[i * 3] = thisPoint[0];
      tetin.pointlist[i * 3 + 1] = thisPoint[1];
      tetin.pointlist[i * 3 + 2] = thisPoint[2];
    }

    // populate input surface mesh
    tetgenio::facet *tetfacet;
    tetgenio::polygon *tetpoly;

    tetin.numberoffacets = breastMeshPoly->GetNumberOfPolys();
    tetin.facetlist = new tetgenio::facet[tetin.numberoffacets];
    //tetin.facetmarkerlist = new int[tetin.numberoffacets];

    for (int i = 0; i < breastMeshPoly->GetNumberOfPolys(); i++)
    {

      vtkCell *thiscell = breastMeshPoly->GetCell(i);

      tetfacet = &tetin.facetlist[i];
      tetfacet->numberofpolygons = 1;
      tetfacet->polygonlist = new tetgenio::polygon[tetfacet->numberofpolygons];
      tetfacet->numberofholes = 0;
      tetfacet->holelist = nullptr;
      tetpoly = &tetfacet->polygonlist[0];
      tetpoly->numberofvertices = thiscell->GetNumberOfPoints();
      tetpoly->vertexlist = new int[tetpoly->numberofvertices];
      for (int j = 0; j < tetpoly->numberofvertices; j++)
      {
        tetpoly->vertexlist[j] = thiscell->GetPointId(j);
        double thisPt[3];
        breastMeshPoly->GetPoint(tetpoly->vertexlist[j], thisPt);
      }
    }

    tetrahedralize(&tetPar, &tetin, &tetout);

    /***************************
     * setup FEBio simulation
     **************************/

    // calculate bounds
    double phantomBounds[6];
    int phantomExtent[6];

    input->GetExtent(phantomExtent);

    if (meshCount == 0)
    {
      // find phantom bounds - only needed on first iteration
      //std::cout << "Searching for breast bounds\n";

#pragma omp sections
      {
#pragma omp section
        // back
        {
          bool foundBreast = false;
          for (int i = phantomExtent[0]; i <= phantomExtent[1]; i++)
          {
            if (!foundBreast)
            {
              for (int j = phantomExtent[2]; j <= phantomExtent[3]; j++)
              {
                for (int k = phantomExtent[4]; k <= phantomExtent[5]; k++)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
                  if (p[0] != tissue.bg && !foundBreast)
                  {
                    double pos[3];
                    int ijk[3] = {i, j, k};
                    input->GetPoint(input->ComputePointId(ijk), pos);
                    phantomBounds[0] = pos[0];
                    foundBreast = true;
                  }
                }
              }
            }
          }
        }

#pragma omp section
        // front
        {
          bool foundBreast = false;
          for (int i = phantomExtent[1]; i >= phantomExtent[0]; i--)
          {
            if (!foundBreast)
            {
              for (int j = phantomExtent[2]; j <= phantomExtent[3]; j++)
              {
                for (int k = phantomExtent[4]; k <= phantomExtent[5]; k++)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
                  if (p[0] != tissue.bg && !foundBreast)
                  {
                    double pos[3];
                    int ijk[3] = {i, j, k};
                    input->GetPoint(input->ComputePointId(ijk), pos);
                    phantomBounds[1] = pos[0];
                    foundBreast = true;
                  }
                }
              }
            }
          }
        }

#pragma omp section
        // left
        {
          bool foundBreast = false;
          for (int j = phantomExtent[2]; j <= phantomExtent[3]; j++)
          {
            if (!foundBreast)
            {
              for (int i = phantomExtent[0]; i <= phantomExtent[1]; i++)
              {
                for (int k = phantomExtent[4]; k <= phantomExtent[5]; k++)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
                  if (p[0] != tissue.bg && !foundBreast)
                  {
                    double pos[3];
                    int ijk[3] = {i, j, k};
                    input->GetPoint(input->ComputePointId(ijk), pos);
                    phantomBounds[2] = pos[1];
                    foundBreast = true;
                  }
                }
              }
            }
          }
        }

#pragma omp section
        // right
        {
          bool foundBreast = false;
          for (int j = phantomExtent[3]; j >= phantomExtent[2]; j--)
          {
            if (!foundBreast)
            {
              for (int i = phantomExtent[0]; i <= phantomExtent[1]; i++)
              {
                for (int k = phantomExtent[4]; k <= phantomExtent[5]; k++)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
                  if (p[0] != tissue.bg && !foundBreast)
                  {
                    double pos[3];
                    int ijk[3] = {i, j, k};
                    input->GetPoint(input->ComputePointId(ijk), pos);
                    phantomBounds[3] = pos[1];
                    foundBreast = true;
                  }
                }
              }
            }
          }
        }

#pragma omp section
        // bottom
        {
          bool foundBreast = false;
          for (int k = phantomExtent[4]; k <= phantomExtent[5]; k++)
          {
            if (!foundBreast)
            {
              for (int i = phantomExtent[0]; i <= phantomExtent[1]; i++)
              {
                for (int j = phantomExtent[2]; j <= phantomExtent[3]; j++)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
                  if (p[0] != tissue.bg && !foundBreast)
                  {
                    double pos[3];
                    int ijk[3] = {i, j, k};
                    input->GetPoint(input->ComputePointId(ijk), pos);
                    phantomBounds[4] = pos[2];
                    foundBreast = true;
                  }
                }
              }
            }
          }
        }

#pragma omp section
        // top
        {
          bool foundBreast = false;
          for (int k = phantomExtent[5]; k >= phantomExtent[4]; k--)
          {
            if (!foundBreast)
            {
              for (int i = phantomExtent[0]; i <= phantomExtent[1]; i++)
              {
                for (int j = phantomExtent[2]; j <= phantomExtent[3]; j++)
                {
                  unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
                  if (p[0] != tissue.bg && !foundBreast)
                  {
                    double pos[3];
                    int ijk[3] = {i, j, k};
                    input->GetPoint(input->ComputePointId(ijk), pos);
                    phantomBounds[5] = pos[2];
                    foundBreast = true;
                  }
                }
              }
            }
          }
        }
      } // end parallel sections
    }

    //Update nipple position each iteration
    //std::cout<<"Searching for nipple position.\n";

    double nipplePos[3];

    bool foundNipple = false;
    for (int i = phantomExtent[1]; i >= phantomExtent[0]; i--)
    {
      if (!foundNipple)
      {
        int voxCount = 0;
        double nippleSum[3] = {0.0, 0.0, 0.0};
        for (int j = phantomExtent[2]; j <= phantomExtent[3]; j++)
        {
          for (int k = phantomExtent[4]; k <= phantomExtent[5]; k++)
          {
            unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(i, j, k));
            // position is mean position in plane farthest from chest wall
            if (p[0] == tissue.nipple)
            {
              double pos[3];
              int ijk[3] = {i, j, k};
              input->GetPoint(input->ComputePointId(ijk), pos);
              voxCount += 1;
              for (int m = 0; m < 3; m++)
              {
                nippleSum[m] += pos[m];
              }
            }
          }
        }
        if (voxCount > 0)
        {
          foundNipple = true;
          for (int m = 0; m < 3; m++)
          {
            nipplePos[m] = nippleSum[m] / voxCount;
          }
        }
      }
    }

    double topPaddleCenter[3];
    double bottomPaddleCenter[3];
    double paddleHeight = 1.0;

    if (meshCount == 0)
    {
      origPaddleWidth = (phantomBounds[3] - phantomBounds[2]) * 1.75;
      origPaddleLength = (phantomBounds[1] - phantomBounds[0]) * 1.75;
    }

    double paddleWidth = origPaddleWidth;
    double paddleLength = origPaddleLength;
    double paddleEdgeRadius = 3.0;
    double paddleLipHeight = 15.0; // must be greater than paddleEdgeRadius

    // paddle retreat after remeshing
    double retreatDist = 10.0 * paddleHeight;

    if (meshCount == 0)
    {
      topPaddleCenter[0] = paddleLength / 2.0 + (phantomBounds[1] - phantomBounds[0]) * 0.125;
      topPaddleCenter[1] = nipplePos[1];
      if (!rotate)
      {
        topPaddleCenter[2] = phantomBounds[5] + 4 * paddleHeight / 2.0;
      }
      else
      {
        topPaddleCenter[2] = phantomBounds[5] + 30 * paddleHeight / 2.0;
      }

      bottomPaddleCenter[0] = paddleLength / 2.0 + (phantomBounds[1] - phantomBounds[0]) * 0.1;
      bottomPaddleCenter[1] = nipplePos[1];
      if (!rotate)
      {
        bottomPaddleCenter[2] = phantomBounds[4] - 4 * paddleHeight / 2.0;
      }
      else
      {
        bottomPaddleCenter[2] = phantomBounds[4] - 30 * paddleHeight / 2.0;
      }
    }
    else
    {
      // use previous interation's paddle position plus a retreat
      for (int i = 0; i < 3; i++)
      {
        topPaddleCenter[i] = prevTopPaddleCenterDerotated[i];
        bottomPaddleCenter[i] = prevBotPaddleCenterDerotated[i];
      }
      topPaddleCenter[2] += retreatDist;
      bottomPaddleCenter[2] -= retreatDist;
    }

    if (autoRemesh)
    {
      // try to complete full compression
      paddleDist = ((topPaddleCenter[2] - bottomPaddleCenter[2] -
                     paddleHeight - thickness) /
                    2.0);
    }
    else
    {
      if (numMeshes >= 3)
      {
        // paddle distance must be greater for every iteration past the first one to
        // account for the distancing done when resuming from the last iteration.
        if (meshCount == 0)
        {
          paddleDist = ((topPaddleCenter[2] - bottomPaddleCenter[2] -
                         paddleHeight - thickness) /
                        2.0) *
                       0.4;
        }
        else if (meshCount == 1)
        {
          paddleDist = paddleDist / 2;
          paddleDist += retreatDist;
        }
        else
        {
          paddleDist -= retreatDist;
          paddleDist = (paddleDist * 2) / (numMeshes - 2);
          paddleDist += retreatDist;
        }
      }
      else
      {
        // 1 or 2 meshes
        if (meshCount == 0)
        {
          paddleDist = ((topPaddleCenter[2] - bottomPaddleCenter[2] -
                         paddleHeight - thickness) /
                        2.0) /
                       numMeshes;
        }
        else if (meshCount == 1)
        {
          paddleDist += retreatDist;
        }
      }
    }

    if (paddleDist < 0.0)
    {
      cerr << "Paddle travel distance error.\n";
      return (1);
    }

    vtkSmartPointer<vtkDoubleArray> bottomPaddleNodes =
        vtkSmartPointer<vtkDoubleArray>::New();

    bottomPaddleNodes->SetNumberOfComponents(3);

    vtkSmartPointer<vtkIntArray> bottomPaddleHexElements =
        vtkSmartPointer<vtkIntArray>::New();

    bottomPaddleHexElements->SetNumberOfComponents(8);

    vtkSmartPointer<vtkIntArray> bottomPaddlePentElements =
        vtkSmartPointer<vtkIntArray>::New();

    bottomPaddlePentElements->SetNumberOfComponents(6);

    vtkSmartPointer<vtkIntArray> bottomPaddleContactSurface =
        vtkSmartPointer<vtkIntArray>::New();

    bottomPaddleContactSurface->SetNumberOfComponents(4);

    vtkSmartPointer<vtkIntArray> breastBackNodes =
        vtkSmartPointer<vtkIntArray>::New();

    breastBackNodes->SetNumberOfComponents(1);

    double paddleDy = paddleWidth / paddleNw;
    double paddleDx = (paddleLength - 2.0 * paddleEdgeRadius) / (paddleNl - 1);
    vtkIdType nodes[8];
    vtkIdType pnodes[6];

    // create main part of bottom paddle
    for (int j = 0; j < paddleNl - 1; j++)
    {
      double xpos = bottomPaddleCenter[0] + paddleLength / 2.0 - j * paddleDx;
      for (int i = 0; i < paddleNw; i++)
      {
        double ypos = bottomPaddleCenter[1] - paddleWidth / 2.0 + i * paddleDy;
        double loc[3];
        // add nodes
        if (i == 0 && j == 0)
        {
          // add all 8 nodes
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[0] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[1] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[3] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[4] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[5] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[7] = bottomPaddleNodes->InsertNextTuple(loc);
        }
        else if (j == 0)
        {
          // add 4 nodes
          nodes[0] = nodes[1];
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[1] = bottomPaddleNodes->InsertNextTuple(loc);
          nodes[3] = nodes[2];
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
          nodes[4] = nodes[5];
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[5] = bottomPaddleNodes->InsertNextTuple(loc);
          nodes[7] = nodes[6];
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = bottomPaddleNodes->InsertNextTuple(loc);
        }
        else if (i == 0)
        {
          // add 4 nodes
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[3] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[7] = bottomPaddleNodes->InsertNextTuple(loc);
          // get left neighbor nodes
          double leftN[8];
          bottomPaddleHexElements->GetTuple((j - 1) * (paddleNw), leftN);
          nodes[0] = (int)leftN[3] - 1;
          nodes[1] = (int)leftN[2] - 1;
          nodes[4] = (int)leftN[7] - 1;
          nodes[5] = (int)leftN[6] - 1;
        }
        else
        {
          // add 2 nodes
          nodes[0] = nodes[1];
          nodes[3] = nodes[2];
          nodes[4] = nodes[5];
          nodes[7] = nodes[6];
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = bottomPaddleNodes->InsertNextTuple(loc);
          // get left neighbor nodes
          double *leftN;
          leftN = bottomPaddleHexElements->GetTuple((j - 1) * (paddleNw) + i);
          nodes[1] = (int)leftN[2] - 1;
          nodes[5] = (int)leftN[6] - 1;
        }

        // add element
        double nodesD[8];
        for (int a = 0; a < 8; a++)
        {
          nodesD[a] = (double)(nodes[a] + 1);
        }
        bottomPaddleHexElements->InsertNextTuple(nodesD);

        // add surface
        bottomPaddleContactSurface->InsertNextTuple4(nodesD[4], nodesD[5], nodesD[6], nodesD[7]);
      }
    }

    int saveTop;

    // curved section
    for (int i = 0; i < paddleNw; i++)
    {
      double ypos = bottomPaddleCenter[1] - paddleWidth / 2.0 + i * paddleDy;
      double xpos = bottomPaddleCenter[0] + paddleLength / 2.0 - (paddleNl - 1) * paddleDx;

      // middle pentahedral
      double leftN[8];
      double loc[3];
      bottomPaddleHexElements->GetTuple((paddleNl - 2) * (paddleNw) + i, leftN);

      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleEdgeRadius;
        pnodes[5] = bottomPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        pnodes[5] = pnodes[0];
      }
      loc[0] = xpos - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleEdgeRadius;
      pnodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
      pnodes[3] = (int)leftN[3] - 1;
      pnodes[4] = (int)leftN[7] - 1;
      pnodes[0] = (int)leftN[2] - 1;
      pnodes[1] = (int)leftN[6] - 1;

      // save 2 nodes for top pent
      int nodSave[2];
      nodSave[0] = pnodes[1];
      nodSave[1] = pnodes[4];

      double pnodesD[6];
      for (int a = 0; a < 6; a++)
      {
        pnodesD[a] = (double)(pnodes[a] + 1);
      }

      bottomPaddlePentElements->InsertNextTuple(pnodesD);

      // bottom pentahedral
      // 0 and 3 the same
      pnodes[1] = pnodes[2];
      pnodes[4] = pnodes[5];
      loc[0] = xpos - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
      pnodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
        pnodes[5] = bottomPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        // from previous edge hex
        pnodes[5] = nodes[1];
      }

      for (int a = 0; a < 6; a++)
      {
        pnodesD[a] = (double)(pnodes[a] + 1);
      }
      bottomPaddlePentElements->InsertNextTuple(pnodesD);

      // bottom edge hex
      nodes[0] = pnodes[5];
      nodes[4] = pnodes[4];
      nodes[1] = pnodes[2];
      nodes[5] = pnodes[1];
      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
        nodes[3] = bottomPaddleNodes->InsertNextTuple(loc);
        loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleEdgeRadius;
        nodes[7] = bottomPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        nodes[3] = nodes[2];
        nodes[7] = nodes[6];
      }

      loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
      nodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
      loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleEdgeRadius;
      nodes[6] = bottomPaddleNodes->InsertNextTuple(loc);

      double nodesD[8];
      for (int a = 0; a < 8; a++)
      {
        nodesD[a] = (double)(nodes[a] + 1);
      }
      bottomPaddleHexElements->InsertNextTuple(nodesD);

      // the edge
      bottomPaddleContactSurface->InsertNextTuple4(nodesD[3], nodesD[7], nodesD[6], nodesD[2]);
      // add the bottom edge as well
      bottomPaddleContactSurface->InsertNextTuple4(nodesD[3], nodesD[2], nodesD[1], nodesD[0]);

      // top pentahedral
      // pnodes 1,4 same
      pnodes[2] = nodSave[0];
      pnodes[5] = nodSave[1];

      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
        pnodes[3] = bottomPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        pnodes[3] = saveTop;
      }

      loc[0] = xpos - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
      pnodes[0] = bottomPaddleNodes->InsertNextTuple(loc);
      saveTop = pnodes[0];

      for (int a = 0; a < 6; a++)
      {
        pnodesD[a] = (double)(pnodes[a] + 1);
      }
      bottomPaddlePentElements->InsertNextTuple(pnodesD);

      // the edge
      bottomPaddleContactSurface->InsertNextTuple4(pnodesD[5], pnodesD[2], pnodesD[0], pnodesD[3]);

      // the fan
      for (int j = 0; j < paddleNh; j++)
      {
        if (j == 0)
        {
          // from top pent
          pnodes[1] = pnodes[0];
          pnodes[4] = pnodes[3];
          pnodes[0] = nodes[5];
          pnodes[3] = nodes[4];
        }
        else
        {
          pnodes[1] = pnodes[2];
          pnodes[4] = pnodes[5];
        }

        double theta = 3.141592654 / 2.0 / (paddleNh) * (j + 1);
        if (j == paddleNh - 1)
        {
          pnodes[2] = nodes[6];
        }
        else
        {
          loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius * sin(theta);
          loc[1] = ypos + paddleDy;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleEdgeRadius + paddleEdgeRadius * cos(theta);
          pnodes[2] = bottomPaddleNodes->InsertNextTuple(loc);
        }
        if (i == 0 && j < paddleNh - 1)
        {
          loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius * sin(theta);
          loc[1] = ypos;
          loc[2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleEdgeRadius + paddleEdgeRadius * cos(theta);
          pnodes[5] = bottomPaddleNodes->InsertNextTuple(loc);
        }
        else if (j < paddleNh - 1)
        {
          double *leftP;

          leftP = bottomPaddlePentElements->GetTuple((i - 1) * (3 + paddleNh) + 3 + j);
          pnodes[5] = (int)leftP[2] - 1;
        }
        else
        {
          pnodes[5] = nodes[7];
        }

        for (int a = 0; a < 6; a++)
        {
          pnodesD[a] = (double)(pnodes[a] + 1);
        }
        bottomPaddlePentElements->InsertNextTuple(pnodesD);

        // the edge
        bottomPaddleContactSurface->InsertNextTuple4(pnodesD[5], pnodesD[4], pnodesD[1], pnodesD[2]);
      }
    }

    // create top paddle nodes
    vtkSmartPointer<vtkDoubleArray> topPaddleNodes =
        vtkSmartPointer<vtkDoubleArray>::New();

    topPaddleNodes->SetNumberOfComponents(3);

    vtkSmartPointer<vtkIntArray> topPaddleHexElements =
        vtkSmartPointer<vtkIntArray>::New();

    topPaddleHexElements->SetNumberOfComponents(8);

    vtkSmartPointer<vtkIntArray> topPaddlePentElements =
        vtkSmartPointer<vtkIntArray>::New();

    topPaddlePentElements->SetNumberOfComponents(6);

    vtkSmartPointer<vtkIntArray> topPaddleContactSurface =
        vtkSmartPointer<vtkIntArray>::New();

    topPaddleContactSurface->SetNumberOfComponents(4);

    // create main part of top paddle
    for (int j = 0; j < paddleNl - 1; j++)
    {
      double xpos = topPaddleCenter[0] + paddleLength / 2.0 - j * paddleDx;
      for (int i = 0; i < paddleNw; i++)
      {
        double ypos = topPaddleCenter[1] - paddleWidth / 2.0 + i * paddleDy;
        double loc[3];
        // add nodes
        if (i == 0 && j == 0)
        {
          // add all 8 nodes
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[0] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[1] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[3] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[4] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[5] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[7] = topPaddleNodes->InsertNextTuple(loc);
        }
        else if (j == 0)
        {
          // add 4 nodes
          nodes[0] = nodes[1];
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[1] = topPaddleNodes->InsertNextTuple(loc);
          nodes[3] = nodes[2];
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = topPaddleNodes->InsertNextTuple(loc);
          nodes[4] = nodes[5];
          loc[0] = xpos;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[5] = topPaddleNodes->InsertNextTuple(loc);
          nodes[7] = nodes[6];
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = topPaddleNodes->InsertNextTuple(loc);
        }
        else if (i == 0)
        {
          // add 4 nodes
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[3] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[7] = topPaddleNodes->InsertNextTuple(loc);
          // get left neighbor nodes
          double leftN[8];
          topPaddleHexElements->GetTuple((j - 1) * (paddleNw), leftN);
          nodes[0] = (int)leftN[3] - 1;
          nodes[1] = (int)leftN[2] - 1;
          nodes[4] = (int)leftN[7] - 1;
          nodes[5] = (int)leftN[6] - 1;
        }
        else
        {
          // add 2 nodes
          nodes[0] = nodes[1];
          nodes[3] = nodes[2];
          nodes[4] = nodes[5];
          nodes[7] = nodes[6];
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
          nodes[2] = topPaddleNodes->InsertNextTuple(loc);
          loc[0] = xpos - paddleDx;
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] + paddleHeight / 2.0;
          nodes[6] = topPaddleNodes->InsertNextTuple(loc);
          // get left neighbor nodes
          double *leftN;
          leftN = topPaddleHexElements->GetTuple((j - 1) * (paddleNw) + i);
          nodes[1] = (int)leftN[2] - 1;
          nodes[5] = (int)leftN[6] - 1;
        }

        // add element
        double nodesD[8];
        for (int a = 0; a < 8; a++)
        {
          nodesD[a] = (double)(nodes[a] + 1);
        }
        topPaddleHexElements->InsertNextTuple(nodesD);

        // add surface
        topPaddleContactSurface->InsertNextTuple4(nodesD[0], nodesD[3], nodesD[2], nodesD[1]);
      }
    }

    int saveBottom;

    // curved section
    for (int i = 0; i < paddleNw; i++)
    {
      double ypos = topPaddleCenter[1] - paddleWidth / 2.0 + i * paddleDy;
      double xpos = topPaddleCenter[0] + paddleLength / 2.0 - (paddleNl - 1) * paddleDx;

      // middle pentahedral
      double leftN[8];
      double loc[3];
      topPaddleHexElements->GetTuple((paddleNl - 2) * (paddleNw) + i, leftN);

      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleEdgeRadius;
        pnodes[5] = topPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        pnodes[5] = pnodes[0];
      }

      loc[0] = xpos - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleEdgeRadius;
      pnodes[2] = topPaddleNodes->InsertNextTuple(loc);
      pnodes[0] = (int)leftN[2] - 1;
      pnodes[1] = (int)leftN[6] - 1;
      pnodes[3] = (int)leftN[3] - 1;
      pnodes[4] = (int)leftN[7] - 1;

      // save 2 nodes for bottom pent
      int nodSave[2];
      nodSave[0] = pnodes[0];
      nodSave[1] = pnodes[3];

      double pnodesD[6];
      for (int a = 0; a < 6; a++)
      {
        pnodesD[a] = (double)(pnodes[a] + 1);
      }

      topPaddlePentElements->InsertNextTuple(pnodesD);

      // top pentahedral
      int tempNode;
      tempNode = pnodes[2];
      pnodes[2] = pnodes[1];
      pnodes[1] = tempNode;
      tempNode = pnodes[5];
      pnodes[5] = pnodes[4];
      pnodes[4] = tempNode;

      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
        pnodes[3] = topPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        pnodes[3] = saveTop;
      }

      loc[0] = xpos - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
      pnodes[0] = topPaddleNodes->InsertNextTuple(loc);
      saveTop = pnodes[0];

      for (int a = 0; a < 6; a++)
      {
        pnodesD[a] = (double)(pnodes[a] + 1);
      }
      topPaddlePentElements->InsertNextTuple(pnodesD);

      // top edge hex
      nodes[0] = pnodes[4];
      nodes[1] = pnodes[1];
      nodes[4] = pnodes[3];
      nodes[5] = pnodes[0];

      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleEdgeRadius;
        nodes[3] = topPaddleNodes->InsertNextTuple(loc);
        loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
        nodes[7] = topPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        nodes[3] = nodes[2];
        nodes[7] = nodes[6];
      }

      loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleEdgeRadius;
      nodes[2] = topPaddleNodes->InsertNextTuple(loc);
      loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
      nodes[6] = topPaddleNodes->InsertNextTuple(loc);

      double nodesD[8];
      for (int a = 0; a < 8; a++)
      {
        nodesD[a] = (double)(nodes[a] + 1);
      }
      topPaddleHexElements->InsertNextTuple(nodesD);

      // the edge
      topPaddleContactSurface->InsertNextTuple4(nodesD[3], nodesD[7], nodesD[6], nodesD[2]);
      // add the bottom edge as well
      topPaddleContactSurface->InsertNextTuple4(nodesD[7], nodesD[4], nodesD[5], nodesD[6]);

      // bottom pentahedral
      // 1 and 4 the same as top pent
      pnodes[0] = nodSave[0];
      pnodes[3] = nodSave[1];

      if (i == 0)
      {
        loc[0] = xpos - paddleEdgeRadius;
        loc[1] = ypos;
        loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
        pnodes[5] = topPaddleNodes->InsertNextTuple(loc);
      }
      else
      {
        // from previous edge hex
        pnodes[5] = saveBottom;
      }

      loc[0] = xpos - paddleEdgeRadius;
      loc[1] = ypos + paddleDy;
      loc[2] = topPaddleCenter[2] - paddleHeight / 2.0;
      pnodes[2] = topPaddleNodes->InsertNextTuple(loc);
      saveBottom = pnodes[2];

      for (int a = 0; a < 6; a++)
      {
        pnodesD[a] = (double)(pnodes[a] + 1);
      }
      topPaddlePentElements->InsertNextTuple(pnodesD);

      // the edge
      topPaddleContactSurface->InsertNextTuple4(pnodesD[5], pnodesD[2], pnodesD[0], pnodesD[3]);

      // the fan
      for (int j = 0; j < paddleNh; j++)
      {
        if (j == 0)
        {
          // from bottom pent
          pnodes[0] = pnodes[1];
          pnodes[3] = pnodes[4];
          // pnodes 2,5 the same
        }
        else
        {
          // pnodes 0,3 stay same
          pnodes[2] = pnodes[1];
          pnodes[5] = pnodes[4];
        }

        // need to set 1 and 4
        double theta = pi / 2.0 / (paddleNh) * (j + 1);
        if (j == paddleNh - 1)
        {
          pnodes[1] = nodes[2];
        }
        else
        {
          loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius * sin(theta);
          loc[1] = ypos + paddleDy;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleEdgeRadius - paddleEdgeRadius * cos(theta);
          pnodes[1] = topPaddleNodes->InsertNextTuple(loc);
        }
        if (i == 0 && j < paddleNh - 1)
        {
          loc[0] = xpos - paddleEdgeRadius - paddleEdgeRadius * sin(theta);
          loc[1] = ypos;
          loc[2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleEdgeRadius - paddleEdgeRadius * cos(theta);
          pnodes[4] = topPaddleNodes->InsertNextTuple(loc);
        }
        else if (j < paddleNh - 1)
        {
          double *leftP;

          leftP = topPaddlePentElements->GetTuple((i - 1) * (3 + paddleNh) + 3 + j);
          pnodes[4] = (int)leftP[1] - 1;
        }
        else
        {
          pnodes[4] = nodes[3];
        }

        for (int a = 0; a < 6; a++)
        {
          pnodesD[a] = (double)(pnodes[a] + 1);
        }
        topPaddlePentElements->InsertNextTuple(pnodesD);

        // the edge
        topPaddleContactSurface->InsertNextTuple4(pnodesD[5], pnodesD[4], pnodesD[1], pnodesD[2]);
      }
    }

    // Initial value before rotation is unchanged.
    double topPaddleCenterRot[3], bottomPaddleCenterRot[3];
    for (int i = 0; i < 3; i++)
    {
      topPaddleCenterRot[i] = topPaddleCenter[i];
      bottomPaddleCenterRot[i] = bottomPaddleCenter[i];
    }

    // rotate paddles
    if (rotate)
    {
      double cRot = cos(angle);
      double sRot = sin(angle);

      for (int i = 0; i < bottomPaddleNodes->GetNumberOfTuples(); i++)
      {
        double t[3];
        double tnew[3];
        bottomPaddleNodes->GetTuple(i, t);
        tnew[0] = t[0];
        tnew[1] = cRot * t[1] - sRot * t[2];
        tnew[2] = sRot * t[1] + cRot * t[2];
        bottomPaddleNodes->SetTuple(i, tnew);
      }

      for (int i = 0; i < topPaddleNodes->GetNumberOfTuples(); i++)
      {
        double t[3];
        double tnew[3];
        topPaddleNodes->GetTuple(i, t);
        tnew[0] = t[0];
        tnew[1] = cRot * t[1] - sRot * t[2];
        tnew[2] = sRot * t[1] + cRot * t[2];
        topPaddleNodes->SetTuple(i, tnew);
      }

      topPaddleCenterRot[0] = topPaddleCenter[0];
      topPaddleCenterRot[1] = cRot * topPaddleCenter[1] - sRot * topPaddleCenter[2];
      topPaddleCenterRot[2] = sRot * topPaddleCenter[1] + cRot * topPaddleCenter[2];

      bottomPaddleCenterRot[0] = bottomPaddleCenter[0];
      bottomPaddleCenterRot[1] = cRot * bottomPaddleCenter[1] - sRot * bottomPaddleCenter[2];
      bottomPaddleCenterRot[2] = sRot * bottomPaddleCenter[1] + cRot * bottomPaddleCenter[2];
    }

    // write FEBio simulation file
    cout << "Writing FEBio simulation file...\n";
    ofstream febFile;
    char febioFilename[256];
    sprintf(febioFilename, "%s/febio_%d_%d.feb", workDir.c_str(), seed, meshCount);
    febFile.open(febioFilename);

    // header
    febFile << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    febFile << "<febio_spec version=\"2.0\">\n";
    febFile << "\t<Module type=\"solid\"/>\n";

    // control
    int timeSteps;
    double dtMin;
    if (autoRemesh)
    {
      timeSteps = 20;
      dtMin = 0.005;
    }
    else
    {
      timeSteps = (int)(ceil(100 / numMeshes));
      dtMin = 0.001;
    }
    double stepSize = 1.0 / ((double)timeSteps);

    febFile << "\t<Control>\n";
    febFile << "\t\t<time_steps>" << timeSteps << "</time_steps>\n";
    febFile << "\t\t<step_size>" << stepSize << "</step_size>\n";
    febFile << "\t\t<max_refs>10</max_refs>\n";
    febFile << "\t\t<max_ups>10</max_ups>\n";
    febFile << "\t\t<dtol>0.001</dtol>\n";
    febFile << "\t\t<etol>0.01</etol>\n";
    febFile << "\t\t<rtol>0</rtol>\n";
    febFile << "\t\t<lstol>0.9</lstol>\n";
    febFile << "\t\t<time_stepper>\n";
    febFile << "\t\t\t<dtmin>" << dtMin << "</dtmin>\n";
    febFile << "\t\t\t<dtmax>" << stepSize << "</dtmax>\n";
    // can lower max_retries to either pass or fail faster.
    if (autoRemesh)
    {
      febFile << "\t\t\t<max_retries>7</max_retries>\n";
    }
    else
    {
      febFile << "\t\t\t<max_retries>50</max_retries>\n";
    }
    febFile << "\t\t\t<opt_iter>10</opt_iter>\n";
    febFile << "\t\t</time_stepper>\n";
    febFile << "\t\t<print_level>PRINT_MAJOR_ITRS</print_level>\n";
    febFile << "\t\t<plot_level>PLOT_MAJOR_ITRS</plot_level>\n";
    febFile << "\t\t<analysis type=\"static\"/>\n";
    febFile << "\t</Control>\n";

    // global variables
    febFile << "\t<Globals>\n";
    febFile << "\t\t<Constants>\n";
    febFile << "\t\t\t<T>0</T>\n";
    febFile << "\t\t\t<R>0</R>\n";
    febFile << "\t\t\t<Fc>0</Fc>\n";
    febFile << "\t\t</Constants>\n";
    febFile << "\t</Globals>\n";

    // Materials
    febFile << "\t<Material>\n";
    febFile << "\t\t<material id=\"1\" name=\"BreastFat\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << fatDensity << "</density>\n";
    febFile << "\t\t\t<E>" << fatModulus << "</E>\n";
    febFile << "\t\t\t<v>" << fatPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"2\" name=\"BreastGland\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << glandDensity << "</density>\n";
    febFile << "\t\t\t<E>" << glandModulus << "</E>\n";
    febFile << "\t\t\t<v>" << glandPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"3\" name=\"BottomPaddleMaterial\" type=\"rigid body\">\n";
    febFile << "\t\t\t<density>1e-05</density>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"4\" name=\"TopPaddleMaterial\" type=\"rigid body\">\n";
    febFile << "\t\t\t<density>1e-05</density>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"5\" name=\"BreastMuscle\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << fatDensity << "</density>\n";
    febFile << "\t\t\t<E>" << fatModulus << "</E>\n";
    febFile << "\t\t\t<v>" << fatPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t\t<material id=\"6\" name=\"BreastMass\" type=\"neo-Hookean\">\n";
    febFile << "\t\t\t<density>" << massDensity << "</density>\n";
    febFile << "\t\t\t<E>" << massModulus << "</E>\n";
    febFile << "\t\t\t<v>" << massPoisson << "</v>\n";
    febFile << "\t\t</material>\n";

    febFile << "\t</Material>\n";

    // Geometry
    febFile << "\t<Geometry>\n";

    // nodes
    febFile << "\t\t<Nodes>\n";

    int numBottomPaddleNodes = bottomPaddleNodes->GetNumberOfTuples();

    // bottom paddle nodes
    for (int i = 0; i < numBottomPaddleNodes; i++)
    {
      double loc[3];
      bottomPaddleNodes->GetTuple(i, loc);
      febFile << "\t\t\t<node id=\"" << i + 1 << "\"> " << loc[0] << ", " << loc[1] << ", " << loc[2] << "</node>\n";
    }

    int numTopPaddleNodes = topPaddleNodes->GetNumberOfTuples();

    // top paddle nodes
    for (int i = 0; i < numTopPaddleNodes; i++)
    {
      double loc[3];
      topPaddleNodes->GetTuple(i, loc);
      febFile << "\t\t\t<node id=\"" << i + 1 + numBottomPaddleNodes << "\"> " << loc[0] << ", " << loc[1] << ", " << loc[2] << "</node>\n";
    }

    // breast nodes
    REAL *breastPt = tetout.pointlist;

    for (int i = 0; i < tetout.numberofpoints; i++)
    {
      febFile << "\t\t\t<node id=\"" << i + 1 + numBottomPaddleNodes + numTopPaddleNodes << "\"> " << *breastPt;
      breastPt++;
      febFile << ", " << *breastPt;
      breastPt++;
      febFile << ", " << *breastPt << "</node>\n";
      breastPt++;
    }

    febFile << "\t\t</Nodes>\n";

    // elements
    int elementCount = 0;

    // bottom paddle hex elements
    febFile << "\t\t<Elements type=\"hex8\" mat=\"3\" elset=\"bottomPaddleHex\">\n";
    for (int i = 0; i < bottomPaddleHexElements->GetNumberOfTuples(); i++)
    {
      double *els;
      els = bottomPaddleHexElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i + 1 + elementCount << "\"> ";
      for (int j = 0; j < 7; j++)
      {
        febFile << (int)(els[j]) << ", ";
      }
      febFile << (int)(els[7]) << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += bottomPaddleHexElements->GetNumberOfTuples();

    // bottom paddle pent elements
    febFile << "\t\t<Elements type=\"penta6\" mat=\"3\" elset=\"bottomPaddlePent\">\n";
    for (int i = 0; i < bottomPaddlePentElements->GetNumberOfTuples(); i++)
    {
      double *els;
      els = bottomPaddlePentElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i + 1 + elementCount << "\"> ";
      for (int j = 0; j < 5; j++)
      {
        febFile << (int)(els[j]) << ", ";
      }
      febFile << (int)(els[5]) << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += bottomPaddlePentElements->GetNumberOfTuples();

    // top paddle hex elements
    febFile << "\t\t<Elements type=\"hex8\" mat=\"4\" elset=\"topPaddleHex\">\n";
    for (int i = 0; i < topPaddleHexElements->GetNumberOfTuples(); i++)
    {
      double *els;
      els = topPaddleHexElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i + 1 + elementCount << "\"> ";
      for (int j = 0; j < 7; j++)
      {
        febFile << (int)(els[j]) + numBottomPaddleNodes << ", ";
      }
      febFile << (int)(els[7]) + numBottomPaddleNodes << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += topPaddleHexElements->GetNumberOfTuples();

    // top paddle pent elements
    febFile << "\t\t<Elements type=\"penta6\" mat=\"4\" elset=\"topPaddlePent\">\n";
    for (int i = 0; i < topPaddlePentElements->GetNumberOfTuples(); i++)
    {
      double *els;
      els = topPaddlePentElements->GetTuple(i);
      febFile << "\t\t\t<elem id=\"" << i + 1 + elementCount << "\"> ";
      for (int j = 0; j < 5; j++)
      {
        febFile << (int)(els[j]) + numBottomPaddleNodes << ", ";
      }
      febFile << (int)(els[5]) + numBottomPaddleNodes << "</elem>\n";
    }
    febFile << "\t\t</Elements>\n";

    elementCount += topPaddlePentElements->GetNumberOfTuples();

    // determine tissue type for breast tets and set muscle node list
    vtkSmartPointer<vtkIntArray> tetTypes =
        vtkSmartPointer<vtkIntArray>::New();

    tetTypes->SetNumberOfComponents(1);
    tetTypes->SetNumberOfTuples(tetout.numberoftetrahedra);

    vtkSmartPointer<vtkIdList> muscleNodes =
        vtkSmartPointer<vtkIdList>::New();

    bool fatFound = false;
    bool muscleFound = false;
    bool glandFound = false;
    bool massFound = false;

    for (int i = 0; i < tetout.numberoftetrahedra; i++)
    {
      double coords[4][3];
      double center[3];

      for (int j = 0; j < 4; j++)
      {
        int idx = tetout.tetrahedronlist[4 * i + j];
        for (int k = 0; k < 3; k++)
        {
          coords[j][k] = tetout.pointlist[idx * 3 + k];
        }
      }

      // calculate center of tet
      vtkTetra::TetraCenter(coords[0], coords[1], coords[2], coords[3], center);

      // find tissue type
      //vtkIdType loc = input->FindPoint(center);
      double pcoords[3];
      int ijk[3];
      int inFOV = input->ComputeStructuredCoordinates(center, ijk, pcoords);

      if (inFOV)
      {
        unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(ijk));

        if (p[0] == tissue.fat || p[0] == tissue.skin)
        {
          // use fat material
          tetTypes->SetTuple1(i, 1);
          fatFound = true;
        }
        else if (p[0] == tissue.muscle)
        {
          // use fat and add nodes to fixed boundary
          tetTypes->SetTuple1(i, 5);
          muscleFound = true;
          for (int j = 0; j < 4; j++)
          {
            muscleNodes->InsertUniqueId(tetout.tetrahedronlist[4 * i + j]);
          }
        }
        else if (p[0] == tissue.mass)
        {
          tetTypes->SetTuple1(i, 6);
          massFound = true;
        }
        else if (p[0] == tissue.gland || p[0] == tissue.TDLU ||
                 p[0] == tissue.duct || p[0] == tissue.nipple || p[0] == tissue.cooper || p[0] == tissue.artery || p[0] == tissue.vein)
        {
          // use gland
          tetTypes->SetTuple1(i, 2);
          glandFound = true;
        }
        else if (p[0] == tissue.bg || p[0] == tissue.paddle)
        {
          tetTypes->SetTuple1(i, 1);
        }
        else
        {
          cerr << "Breast tetrahedron center material = " << p[0] << " unrecognized\n";
          return (1);
        }
      }
      else
      {
        // center outside FOV, make it fat
        tetTypes->SetTuple1(i, 1);
      }
    }

    // breast fat elements
    if (fatFound)
    {
      febFile << "\t\t<Elements type=\"tet4\" mat=\"1\" elset=\"breastFatTet\">\n";
      for (int i = 0; i < tetout.numberoftetrahedra; i++)
      {
        if ((int)tetTypes->GetTuple1(i) == 1)
        {
          febFile << "\t\t\t<elem id=\"" << elementCount + 1 << "\"> ";
          int n;
          for (int j = 0; j < 3; j++)
          {
            febFile << tetout.tetrahedronlist[4 * i + j] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4 * i + 3] + numBottomPaddleNodes + numTopPaddleNodes + 1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // breast gland elements
    if (glandFound)
    {
      febFile << "\t\t<Elements type=\"tet4\" mat=\"2\" elset=\"breastGlandTet\">\n";
      for (int i = 0; i < tetout.numberoftetrahedra; i++)
      {
        if ((int)tetTypes->GetTuple1(i) == 2)
        {
          febFile << "\t\t\t<elem id=\"" << elementCount + 1 << "\"> ";
          int n;
          for (int j = 0; j < 3; j++)
          {
            febFile << tetout.tetrahedronlist[4 * i + j] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4 * i + 3] + numBottomPaddleNodes + numTopPaddleNodes + 1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // mass elements
    if (massFound)
    {
      febFile << "\t\t<Elements type=\"tet4\" mat=\"6\" elset=\"breastMassTet\">\n";
      for (int i = 0; i < tetout.numberoftetrahedra; i++)
      {
        if ((int)tetTypes->GetTuple1(i) == 6)
        {
          febFile << "\t\t\t<elem id=\"" << elementCount + 1 << "\"> ";
          int n;
          for (int j = 0; j < 3; j++)
          {
            febFile << tetout.tetrahedronlist[4 * i + j] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4 * i + 3] + numBottomPaddleNodes + numTopPaddleNodes + 1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // muscle elements
    if (muscleFound)
    {
      febFile << "\t\t<Elements type=\"tet4\" mat=\"5\" elset=\"breastMuscleTet\">\n";
      for (int i = 0; i < tetout.numberoftetrahedra; i++)
      {
        if ((int)tetTypes->GetTuple1(i) == 5)
        {
          febFile << "\t\t\t<elem id=\"" << elementCount + 1 << "\"> ";
          int n;
          for (int j = 0; j < 3; j++)
          {
            febFile << tetout.tetrahedronlist[4 * i + j] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", ";
          }
          febFile << tetout.tetrahedronlist[4 * i + 3] + numBottomPaddleNodes + numTopPaddleNodes + 1 << "</elem>\n";
          elementCount += 1;
        }
      }
      febFile << "\t\t</Elements>\n";
    }

    // end geometry
    febFile << "\t</Geometry>\n";

    // boundary - muscle tissue fixed
    febFile << "\t<Boundary>\n";

    febFile << "\t\t<fix bc=\"xyz\">\n";
    for (int i = 0; i < muscleNodes->GetNumberOfIds(); i++)
    {
      vtkIdType n;
      n = muscleNodes->GetId(i);
      febFile << "\t\t\t<node id=\"" << (n + numBottomPaddleNodes + numTopPaddleNodes + 1) << "\"/>\n";
    }
    febFile << "\t\t</fix>\n";

    febFile << "\t\t<fix bc=\"uvw\">\n";
    for (int i = 0; i < muscleNodes->GetNumberOfIds(); i++)
    {
      vtkIdType n;
      n = muscleNodes->GetId(i);
      febFile << "\t\t\t<node id=\"" << (n + numBottomPaddleNodes + numTopPaddleNodes + 1) << "\"/>\n";
    }
    febFile << "\t\t</fix>\n";

    febFile << "\t</Boundary>\n";

    // Gravity
    febFile << "\t<Loads>\n";
    febFile << "\t\t<body_load type=\"const\">\n";
    febFile << "\t\t\t<z lc=\"1\">4905</z>\n"; //gravity in kg*mm/s^2 assuming 0.5 kg breast
    febFile << "\t\t</body_load>\n";
    febFile << "\t</Loads>\n";

    // Contacts

    // bottom contact
    febFile << "\t<Contact>\n";
    febFile << "\t\t<contact type=\"facet-to-facet sliding\">\n";
    febFile << "\t\t\t<laugon>0</laugon>\n";
    febFile << "\t\t\t<tolerance>1.0</tolerance>\n";
    febFile << "\t\t\t<penalty>25</penalty>\n";
    febFile << "\t\t\t<two_pass>0</two_pass>\n";
    febFile << "\t\t\t<auto_penalty>1</auto_penalty>\n";
    febFile << "\t\t\t<fric_coeff>0</fric_coeff>\n";
    febFile << "\t\t\t<fric_penalty>0</fric_penalty>\n";
    febFile << "\t\t\t<search_tol>0.02</search_tol>\n";
    febFile << "\t\t\t<minaug>0</minaug>\n";
    febFile << "\t\t\t<maxaug>10</maxaug>\n";
    febFile << "\t\t\t<gaptol>0</gaptol>\n";
    febFile << "\t\t\t<seg_up>0</seg_up>\n";
    febFile << "\t\t\t<surface type=\"master\">\n";
    for (int i = 0; i < bottomPaddleContactSurface->GetNumberOfTuples(); i++)
    {
      double *n;
      n = bottomPaddleContactSurface->GetTuple4(i);
      febFile << "\t\t\t\t<quad4 id=\"" << i + 1 << "\"> " << (int)(n[0]) << ", " << (int)(n[1]) << ", " << (int)(n[2]) << ", " << (int)(n[3]) << "</quad4>\n";
    }
    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t\t<surface type=\"slave\">\n";
    int triCount = 1;

    double zThresh = nipplePos[2];

    for (int i = 0; i < tetout.numberoftrifaces; i++)
    {
      // test if triface is bottom contact
      double meanZ = 0.0;
      double meanX = 0.0;
      double meanY = 0.0;
      for (int j = 0; j < 3; j++)
      {
        meanZ += tetout.pointlist[3 * (tetout.trifacelist[3 * i + j]) + 2];
        meanY += tetout.pointlist[3 * (tetout.trifacelist[3 * i + j]) + 1];
        meanX += tetout.pointlist[3 * (tetout.trifacelist[3 * i + j])];
      }
      meanZ = meanZ / 3.0;
      meanY = meanY / 3.0;
      meanX = meanX / 3.0;

      if (rotate)
      {
        meanZ = -sin(angle) * meanY + cos(angle) * meanZ;
        zThresh = -sin(angle) * nipplePos[1] + cos(angle) * nipplePos[2];
      }

      if (meanZ < zThresh)
      {
        febFile << "\t\t\t\t<tri3 id=\"" << triCount << "\"> " << tetout.trifacelist[3 * i] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", " << tetout.trifacelist[3 * i + 1] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", " << tetout.trifacelist[3 * i + 2] + numBottomPaddleNodes + numTopPaddleNodes + 1 << "</tri3>\n";
        triCount += 1;
      }
    }

    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t</contact>\n";

    // top contact
    febFile << "\t\t<contact type=\"facet-to-facet sliding\">\n";
    febFile << "\t\t\t<laugon>0</laugon>\n";
    febFile << "\t\t\t<tolerance>1.0</tolerance>\n";
    febFile << "\t\t\t<penalty>25</penalty>\n";
    febFile << "\t\t\t<two_pass>0</two_pass>\n";
    febFile << "\t\t\t<auto_penalty>1</auto_penalty>\n";
    febFile << "\t\t\t<fric_coeff>0</fric_coeff>\n";
    febFile << "\t\t\t<fric_penalty>0</fric_penalty>\n";
    febFile << "\t\t\t<search_tol>0.02</search_tol>\n";
    febFile << "\t\t\t<minaug>0</minaug>\n";
    febFile << "\t\t\t<maxaug>10</maxaug>\n";
    febFile << "\t\t\t<gaptol>0</gaptol>\n";
    febFile << "\t\t\t<seg_up>0</seg_up>\n";
    febFile << "\t\t\t<surface type=\"master\">\n";
    for (int i = 0; i < topPaddleContactSurface->GetNumberOfTuples(); i++)
    {
      double *n;
      n = topPaddleContactSurface->GetTuple4(i);
      febFile << "\t\t\t\t<quad4 id=\"" << i + 1 << "\"> " << (int)(n[0]) + numBottomPaddleNodes << ", " << (int)(n[1]) + numBottomPaddleNodes << ", " << (int)(n[2]) + numBottomPaddleNodes << ", " << (int)(n[3]) + numBottomPaddleNodes << "</quad4>\n";
    }
    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t\t<surface type=\"slave\">\n";
    triCount = 1;

    for (int i = 0; i < tetout.numberoftrifaces; i++)
    {
      // test if triface is top contact
      double meanZ = 0.0;
      double meanX = 0.0;
      double meanY = 0.0;
      for (int j = 0; j < 3; j++)
      {
        meanZ += tetout.pointlist[3 * (tetout.trifacelist[3 * i + j]) + 2];
        meanY += tetout.pointlist[3 * (tetout.trifacelist[3 * i + j]) + 1];
        meanX += tetout.pointlist[3 * (tetout.trifacelist[3 * i + j])];
      }
      meanZ = meanZ / 3.0;
      meanY = meanY / 3.0;
      meanX = meanX / 3.0;

      if (rotate)
      {
        meanZ = -sin(angle) * meanY + cos(angle) * meanZ;
        zThresh = -sin(angle) * nipplePos[1] + cos(angle) * nipplePos[2];
      }

      if (meanZ >= zThresh)
      {
        febFile << "\t\t\t\t<tri3 id=\"" << triCount << "\"> " << tetout.trifacelist[3 * i] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", " << tetout.trifacelist[3 * i + 1] + numBottomPaddleNodes + numTopPaddleNodes + 1 << ", " << tetout.trifacelist[3 * i + 2] + numBottomPaddleNodes + numTopPaddleNodes + 1 << "</tri3>\n";
        triCount += 1;
      }
    }

    febFile << "\t\t\t</surface>\n";

    febFile << "\t\t</contact>\n";

    febFile << "\t</Contact>\n";

    // constraints
    febFile << "\t<Constraints>\n";
    febFile << "\t\t<rigid_body mat=\"4\">\n";
    febFile << "\t\t\t<fixed bc=\"x\"/>\n";
    if (!rotate)
    {
      febFile << "\t\t\t<fixed bc=\"y\"/>\n";
    }
    febFile << "\t\t\t<fixed bc=\"Rx\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Ry\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Rz\"/>\n";
    febFile << "\t\t</rigid_body>\n";
    febFile << "\t\t<rigid_body mat=\"3\">\n";
    febFile << "\t\t\t<fixed bc=\"x\"/>\n";
    if (!rotate)
    {
      febFile << "\t\t\t<fixed bc=\"y\"/>\n";
    }
    febFile << "\t\t\t<fixed bc=\"Rx\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Ry\"/>\n";
    febFile << "\t\t\t<fixed bc=\"Rz\"/>\n";
    febFile << "\t\t</rigid_body>\n";
    febFile << "\t\t<rigid_body mat=\"4\">\n";
    febFile << "\t\t\t<prescribed bc=\"z\" lc=\"2\">1</prescribed>\n";
    febFile << "\t\t</rigid_body>\n";
    if (rotate)
    {
      febFile << "\t\t<rigid_body mat=\"4\">\n";
      febFile << "\t\t\t<prescribed bc=\"y\" lc=\"4\">1</prescribed>\n";
      febFile << "\t\t</rigid_body>\n";
    }
    febFile << "\t\t<rigid_body mat=\"3\">\n";
    febFile << "\t\t\t<prescribed bc=\"z\" lc=\"3\">1</prescribed>\n";
    febFile << "\t\t</rigid_body>\n";
    if (rotate)
    {
      febFile << "\t\t<rigid_body mat=\"3\">\n";
      febFile << "\t\t\t<prescribed bc=\"y\" lc=\"5\">1</prescribed>\n";
      febFile << "\t\t</rigid_body>\n";
    }
    febFile << "\t</Constraints>\n";

    // load curves
    febFile << "\t<LoadData>\n";
    febFile << "\t\t<loadcurve id=\"1\" type=\"linear\">\n";
    febFile << "\t\t\t<point>0,0</point>\n";

    // currently have kept gravity very low but non-zero for stability.
    //    if(meshCount == 0)
    //    {
    //      febFile << "\t\t\t<point>0.3,1.0</point>\n";
    //      febFile << "\t\t\t<point>1,1.0</point>\n";
    //    }
    //    else
    //    {

    //      febFile << "\t\t\t<point>0.3,0</point>\n";
    //      febFile << "\t\t\t<point>0.4,1.0</point>\n";

    // Left a miniscule force so that the breast wouldn't be underconstrained
    // (causing problem to diverge).
    febFile << "\t\t\t<point>1,0.001</point>\n";
    //    }

    febFile << "\t\t</loadcurve>\n";
    if (rotate)
    {
      febFile << "\t\t<loadcurve id=\"2\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << -paddleDist * cos(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"3\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << paddleDist * cos(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"4\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << paddleDist * sin(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"5\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1," << -paddleDist * sin(angle) << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
    }
    else
    {
      // meet in center, each paddle travels an equal distance
      febFile << "\t\t<loadcurve id=\"2\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1, " << -paddleDist << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
      febFile << "\t\t<loadcurve id=\"3\" type=\"linear\">\n";
      febFile << "\t\t\t<point>0,0</point>\n";
      febFile << "\t\t\t<point>0.1,0</point>\n";
      febFile << "\t\t\t<point>1, " << paddleDist << "</point>\n";
      febFile << "\t\t</loadcurve>\n";
    }
    febFile << "\t</LoadData>\n";

    // output
    febFile << "\t<Output>\n";
    febFile << "\t\t<plotfile type=\"febio\">\n";
    febFile << "\t\t\t<var type=\"displacement\"/>\n";
    febFile << "\t\t\t<compression>0</compression>\n";
    febFile << "\t\t</plotfile>\n";
    febFile << "\t</Output>\n";

    // end feb file
    febFile << "</febio_spec>\n";
    febFile.close();

    /************
     * run FEBio
     ***********/

    int FEBioResult;

    cout << "Running FEBio compression modeling...\n";

    char febioOutFilename[256];
    sprintf(febioOutFilename, "%s/febio_%d_%d.xplt", workDir.c_str(), seed, meshCount);
    char febioCmd[512];
    sprintf(febioCmd, "febio2 -i %s -p %s -nosplash > /dev/null", febioFilename, febioOutFilename);

    FEBioResult = system(febioCmd);

    if (FEBioResult == -1)
    {
      std::cerr << "febioCmd Error\n";
    }

    /**********************
     * Parse FEBio result
     *********************/

    // extract node displacements for final iteration
    cout << "Generating voxelized compressed breast...\n";

    FILE *fefile;
    uint32_t dword;
    long timePos;
    long oldTimePos;
    int seekFail = 0;
    int readGood = 1;

    size_t numRead;

    fefile = fopen(febioOutFilename, "rb");

    if (fefile == NULL)
    {
      cerr << "Cannot read FEBio output file\n";
      return (1);
    }

    // check if file format correct
    numRead = fread(&dword, 1, 4, fefile);
    if (dword != 0x00464542)
    {
      cerr << "FEBio file invalid format\n";
      return (1);
    }

    // root block ID
    numRead = fread(&dword, 1, 4, fefile);
    // root block size
    numRead = fread(&dword, 1, 4, fefile);
    // skip root block
    fseek(fefile, dword, SEEK_CUR);

    // find final time step

    // state id
    numRead = fread(&dword, 1, 4, fefile);
    // state size
    numRead = fread(&dword, 1, 4, fefile);
    // save current time step location
    timePos = ftell(fefile);

    seekFail = fseek(fefile, dword, SEEK_CUR);

    while (!seekFail && readGood)
    {
      oldTimePos = timePos;
      // state id
      readGood = fread(&dword, 1, 4, fefile);
      // state size
      readGood = fread(&dword, 1, 4, fefile);
      // save current time step location
      timePos = ftell(fefile);
      seekFail = fseek(fefile, dword, SEEK_CUR);
    }

    // seek to start of last time step
    fseek(fefile, oldTimePos, SEEK_SET);

    // state header
    numRead = fread(&dword, 1, 4, fefile);

    // state header size
    numRead = fread(&dword, 1, 4, fefile);

    // state time tag
    numRead = fread(&dword, 1, 4, fefile);

    // state time size
    numRead = fread(&dword, 1, 4, fefile);

    float stateTime;

    numRead = fread(&stateTime, 1, dword, fefile);

    // state data tag
    numRead = fread(&dword, 1, 4, fefile);

    // state data size
    numRead = fread(&dword, 1, 4, fefile);

    // state node data tag
    numRead = fread(&dword, 1, 4, fefile);

    // node data size
    numRead = fread(&dword, 1, 4, fefile);

    // state variable tag
    numRead = fread(&dword, 1, 4, fefile);

    // node data size
    numRead = fread(&dword, 1, 4, fefile);

    // region id tag = 1 should be 0 according to spec
    numRead = fread(&dword, 1, 4, fefile);

    // region id size
    numRead = fread(&dword, 1, 4, fefile);

    // region id
    numRead = fread(&dword, 1, 4, fefile);

    // data tag
    numRead = fread(&dword, 1, 4, fefile);

    uint32_t dispDataSize;

    // data size
    numRead = fread(&dword, 1, 4, fefile);

    // 0 - not in xplt spec
    numRead = fread(&dword, 1, 4, fefile);
    // the displacement data size (bytes)
    numRead = fread(&dispDataSize, 1, 4, fefile);

    // allocate memory
    float *rawDisp;

    rawDisp = new float[dispDataSize / 4];
    numRead = fread(rawDisp, sizeof *rawDisp, dispDataSize / 4, fefile);

    if (numRead == dispDataSize / 4)
    {
      fclose(fefile);
    }
    else
    {
      fclose(fefile);
      cerr << "Data read failure\n";
      return (1);
    }

    vtkSmartPointer<vtkFloatArray> displace =
        vtkSmartPointer<vtkFloatArray>::New();

    displace->SetNumberOfComponents(3);
    int numNodes = dispDataSize / 4 / 3;
    displace->SetNumberOfTuples(numNodes);

    float *dispPtr = rawDisp;

    for (int i = 0; i < numNodes; i++)
    {
      displace->SetTuple(i, dispPtr);
      dispPtr += 3;
    }

    double topPos[3];
    double bottomPos[3];
    double topDisp[3];
    double bottomDisp[3];
    for (int i = 0; i < 3; i++)
    {
      bottomDisp[i] = rawDisp[i];
      bottomPos[i] = bottomPaddleCenterRot[i] + bottomDisp[i];
      topDisp[i] = rawDisp[numBottomPaddleNodes * 3 + i];
      topPos[i] = topPaddleCenterRot[i] + topDisp[i];
    }

    // derotate paddles post compression, for use in next iteration.
    prevTopPaddleCenterDerotated[0] = topPos[0];
    prevTopPaddleCenterDerotated[1] = cos(angle) * topPos[1] + sin(angle) * topPos[2];
    prevTopPaddleCenterDerotated[2] = -sin(angle) * topPos[1] + cos(angle) * topPos[2];

    prevBotPaddleCenterDerotated[0] = bottomPos[0];
    prevBotPaddleCenterDerotated[1] = cos(angle) * bottomPos[1] + sin(angle) * bottomPos[2];
    prevBotPaddleCenterDerotated[2] = -sin(angle) * bottomPos[1] + cos(angle) * bottomPos[2];

    double finalThickness;

    // got rid of pow(topPos[0]-bottomPos[0],2.0)+ since don't care about
    // about the x distance, which may not be zero, as we are after the thickness
    // in between the paddles in the y-z plane (shortest distance).
    finalThickness = sqrt(pow(topPos[1] - bottomPos[1], 2.0) +
                          pow(topPos[2] - bottomPos[2], 2.0)) -
                     paddleHeight;

    // how much further we went than last time
    double gain = minimumFinalThicknessSoFar - finalThickness;
    if (gain > 0.0)
    {
      minimumFinalThicknessSoFar = finalThickness;
    }

    cout << "Compressed breast thickness = " << finalThickness << " mm\n";

    if (finalThickness - thickness <= 0.1 || (finalThickness - thickness <= 5.0 && gain >= 0 && gain <= 1))
    {
      done = true;
    }

    delete[] rawDisp;

    // have displacement data

    // find nodes that escaped the input voxel domain
    // expand extent if necessary
    int compressedExtent[6];
    int inputExtent[6];
    input->GetExtent(inputExtent);
    input->GetExtent(compressedExtent);
    input->GetSpacing(inputSpacing);
    input->GetOrigin(inputOrigin);

    double compressedBounds[6];

    for (int i = 0; i < 3; i++)
    {
      inputBounds[2 * i] = inputOrigin[i] + inputExtent[2 * i] * inputSpacing[i] - inputSpacing[i] / 2.0;
      inputBounds[2 * i + 1] = inputOrigin[i] + inputExtent[2 * i + 1] * inputSpacing[i] + inputSpacing[i] / 2.0;
      compressedBounds[2 * i] = inputBounds[2 * i];
      compressedBounds[2 * i + 1] = inputBounds[2 * i + 1];
    }

    for (int i = 0; i < tetout.numberofpoints; i++)
    {
      int myId = i + 1 + numBottomPaddleNodes + numTopPaddleNodes;
      double finalCoords[3];
      double dr[3];

      displace->GetTuple(myId, dr);

      for (int j = 0; j < 3; j++)
      {
        finalCoords[j] = tetout.pointlist[i * 3 + j] + dr[j];
      }

      for (int j = 0; j < 3; j++)
      {
        if (fabs(dr[j]) < 200.0)
        {
          if (finalCoords[j] < compressedBounds[2 * j])
          {
            compressedBounds[2 * j] = finalCoords[j];
          }
          else if (finalCoords[j] > compressedBounds[2 * j + 1])
          {
            compressedBounds[2 * j + 1] = finalCoords[j];
          }
        }
        else
        {
          cout << "id = " << myId << " dr = " << dr[0] << " " << dr[1] << " " << dr[2] << " coords = " << finalCoords[0] << " " << finalCoords[1] << " " << finalCoords[2] << "\n";
        }
      }
    }

    // new extent - add 5 mm buffer if changed
    for (int i = 0; i < 3; i++)
    {
      if (compressedBounds[2 * i] < inputBounds[2 * i])
      {
        compressedExtent[2 * i] = (int)floor((compressedBounds[2 * i] - inputOrigin[i]) / inputSpacing[i]) - (int)ceil(5.0 / inputSpacing[i]);
      }
      else
      {
        compressedExtent[2 * i] = inputExtent[2 * i];
      }
      if (compressedBounds[2 * i + 1] > inputBounds[2 * i + 1])
      {
        compressedExtent[2 * i + 1] = (int)ceil((compressedBounds[2 * i + 1] - inputOrigin[i]) / inputSpacing[i]) + (int)ceil(5.0 / inputSpacing[i]);
      }
      else
      {
        compressedExtent[2 * i + 1] = inputExtent[2 * i + 1];
      }
    }

    //cout << "old bounds = " << inputBounds[0] << " " << inputBounds[1] << " " << inputBounds[2] << " " <<
    //  inputBounds[3] << " " << inputBounds[4] << " " << inputBounds[5] << "\n";

    //cout << "new bounds = " << compressedBounds[0] << " " << compressedBounds[1] << " " << compressedBounds[2] << " " <<
    //  compressedBounds[3] << " " << compressedBounds[4] << " " << compressedBounds[5] << "\n";

    // setup compressed phantom
    vtkSmartPointer<vtkImageData> breastCompressed =
        vtkSmartPointer<vtkImageData>::New();

    breastCompressed->SetSpacing(inputSpacing);

    breastCompressed->SetExtent(compressedExtent);

    //cout << "old extent = " << inputExtent[0] << " " << inputExtent[1] << " " << inputExtent[2] << " " <<
    //  inputExtent[3] << " " << inputExtent[4] << " " << inputExtent[5] << "\n";

    //cout << "new extent = " << compressedExtent[0] << " " << compressedExtent[1] << " " << compressedExtent[2] << " " <<
    //  compressedExtent[3] << " " << compressedExtent[4] << " " << compressedExtent[5] << "\n";

    breastCompressed->SetOrigin(input->GetOrigin());

    breastCompressed->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

    // initialized to background tissue
    unsigned char *voxVal = static_cast<unsigned char *>(breastCompressed->GetScalarPointer());
    //unsigned char* voxValDebug = static_cast<unsigned char *>(breastCompressed->GetScalarPointer());
    int dim[3];
    breastCompressed->GetDimensions(dim);
    const long long int numElements = dim[0] * dim[1] * dim[2];
    for (long long int i = 0; i < numElements; i++)
    {
      voxVal[i] = tissue.bg;
    }

    int breastCompressedExtent[6];
    breastCompressed->GetExtent(breastCompressedExtent);

    double breastCompressedBounds[6];
    breastCompressed->GetBounds(breastCompressedBounds);

    // interate over breast tets, labeling voxels
    double spacing[3];
    breastCompressed->GetSpacing(spacing);
    double origin[3];
    breastCompressed->GetOrigin(origin);

    std::cout << "Assigning Tissues\n";

    vtkSmartPointer<vtkPoints> newLesionLoc = vtkSmartPointer<vtkPoints>::New();

    // translate lesion locations
    for (int i = 0; i < tetout.numberoftetrahedra; i++)
    {

      vtkSmartPointer<vtkTetra> myTet =
          vtkSmartPointer<vtkTetra>::New();

      double *parCoords = myTet->GetParametricCoords();

      double initialCoords[4][3];
      double finalCoords[4][3];
      double displacements[4][3];
      int myNodes[4];

      for (int j = 0; j < 4; j++)
      {
        myNodes[j] = tetout.tetrahedronlist[4 * i + j];
        double dr[3];
        displace->GetTuple(myNodes[j] + numBottomPaddleNodes + numTopPaddleNodes, dr);
        for (int k = 0; k < 3; k++)
        {
          initialCoords[j][k] = tetout.pointlist[myNodes[j] * 3 + k];
          displacements[j][k] = dr[k];
          finalCoords[j][k] = initialCoords[j][k] + displacements[j][k];
          parCoords[4 * j + k] = initialCoords[j][k] + displacements[j][k];
        }
      }

      // check lesion insertion points to see if in tetrahedron - initialCoords
      double dummyPt[3] = {999.0, 999.0, 999.0};

      for (int j = 0; j < numLesionLoc; j++)
      {
        double pt[3];
        lesionLoc->GetPoint(j, pt);

        if (pt[0] != 999.0)
        {
          double bcoords[4];
          vtkTetra::BarycentricCoords(pt, initialCoords[0], initialCoords[1],
                                      initialCoords[2], initialCoords[3], bcoords);
          if (bcoords[0] >= 0.0 && bcoords[1] >= 0.0 && bcoords[2] >= 0.0 && bcoords[3] >= 0.0)
          {

            // inside tet - calculate displaced point
            double disp[3];

            // interpolate displacements
            for (int d = 0; d < 3; d++)
            {
              disp[d] = bcoords[0] * displacements[0][d] + bcoords[1] * displacements[1][d] +
                        bcoords[2] * displacements[2][d] + bcoords[3] * displacements[3][d];
            }

            double newLoc[3] = {pt[0] + disp[0], pt[1] + disp[1], pt[2] + disp[2]};

            newLesionLoc->InsertNextPoint(newLoc);

            // remove point from search list
            lesionLoc->SetPoint(j, dummyPt);
          }
        }
      }
    }

#pragma omp parallel for
    for (int i = 0; i < tetout.numberoftetrahedra; i++)
    {
      vtkSmartPointer<vtkTetra> myTet =
          vtkSmartPointer<vtkTetra>::New();

      double *parCoords = myTet->GetParametricCoords();
      double initialCoords[4][3];
      double finalCoords[4][3];
      double displacements[4][3];
      int myNodes[4];

      for (int j = 0; j < 4; j++)
      {
        myNodes[j] = tetout.tetrahedronlist[4 * i + j];
        double dr[3];
        displace->GetTuple(myNodes[j] + numBottomPaddleNodes + numTopPaddleNodes, dr);
        for (int k = 0; k < 3; k++)
        {
          initialCoords[j][k] = tetout.pointlist[myNodes[j] * 3 + k];
          displacements[j][k] = dr[k];
          finalCoords[j][k] = initialCoords[j][k] + displacements[j][k];
          parCoords[4 * j + k] = initialCoords[j][k] + displacements[j][k];
        }
      }

      // setup bounding box
      int ijk[3];
      double pcoords[3];

      breastCompressed->ComputeStructuredCoordinates(finalCoords[0], ijk, pcoords);
      int bounds[6] = {ijk[0], ijk[0], ijk[1], ijk[1], ijk[2], ijk[2]};

      for (int j = 1; j < 4; j++)
      {
        breastCompressed->ComputeStructuredCoordinates(finalCoords[j], ijk, pcoords);
        if (ijk[0] < bounds[0])
        {
          bounds[0] = ijk[0];
        }
        if (ijk[0] > bounds[1])
        {
          bounds[1] = ijk[0];
        }
        if (ijk[1] < bounds[2])
        {
          bounds[2] = ijk[1];
        }
        if (ijk[1] > bounds[3])
        {
          bounds[3] = ijk[1];
        }
        if (ijk[2] < bounds[4])
        {
          bounds[4] = ijk[2];
        }
        if (ijk[2] > bounds[5])
        {
          bounds[5] = ijk[2];
        }
      }

      for (int m = 0; m < 6; m += 2)
      {
        if (bounds[m] < breastCompressedExtent[m])
        {
          bounds[m] = breastCompressedExtent[m];
        }
        if (bounds[m + 1] > breastCompressedExtent[m + 1])
        {
          bounds[m + 1] = breastCompressedExtent[m + 1];
        }
      }

      // label voxels inside tet
      for (int a = bounds[0]; a <= bounds[1]; a++)
      {
        for (int b = bounds[2]; b <= bounds[3]; b++)
        {
          for (int c = bounds[4]; c <= bounds[5]; c++)
          {
            double bcoords[4];
            double disp[3];
            // inside tet?
            double loc[3] = {origin[0] + a * spacing[0], origin[1] + b * spacing[1],
                             origin[2] + c * spacing[2]};

            vtkTetra::BarycentricCoords(loc, finalCoords[0], finalCoords[1],
                                        finalCoords[2], finalCoords[3], bcoords);

            if (bcoords[0] >= 0.0 && bcoords[1] >= 0.0 && bcoords[2] >= 0.0 && bcoords[3] >= 0.0)
            {

              // inside tet
              // interpolate displacements
              for (int d = 0; d < 3; d++)
              {
                disp[d] = bcoords[0] * displacements[0][d] + bcoords[1] * displacements[1][d] +
                          bcoords[2] * displacements[2][d] + bcoords[3] * displacements[3][d];
              }

              // location in uncompressed phantom
              double oldLoc[3] = {loc[0] - disp[0], loc[1] - disp[1], loc[2] - disp[2]};

              // voxel in uncompressed phantom
              int inOrig = input->ComputeStructuredCoordinates(oldLoc, ijk, pcoords);

              if (inOrig)
              {
                unsigned char *p = static_cast<unsigned char *>(input->GetScalarPointer(ijk));
                unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(a, b, c));
                *q = *p;
              }
            }
          }
        }
      }
    }
    //}
    //omp_set_dynamic(1);

    if (done)
    {
      std::cout << "Labelling/inserting paddles\n";
      // label paddles
      double crot = cos(angle);
      double srot = sin(angle);

      // bounding box for top paddle
      int topBound[6];

      double topCorners[8][3];
      // unrotated

      topCorners[0][0] = topPaddleCenter[0] - paddleLength / 2.0;
      topCorners[1][0] = topPaddleCenter[0] - paddleLength / 2.0;
      topCorners[2][0] = topPaddleCenter[0] - paddleLength / 2.0;
      topCorners[3][0] = topPaddleCenter[0] - paddleLength / 2.0;
      topCorners[4][0] = topPaddleCenter[0] + paddleLength / 2.0;
      topCorners[5][0] = topPaddleCenter[0] + paddleLength / 2.0;
      topCorners[6][0] = topPaddleCenter[0] + paddleLength / 2.0;
      topCorners[7][0] = topPaddleCenter[0] + paddleLength / 2.0;

      topCorners[0][1] = topPaddleCenter[1] - paddleWidth / 2.0;
      topCorners[1][1] = topPaddleCenter[1] + paddleWidth / 2.0;
      topCorners[2][1] = topPaddleCenter[1] + paddleWidth / 2.0;
      topCorners[3][1] = topPaddleCenter[1] - paddleWidth / 2.0;
      topCorners[4][1] = topPaddleCenter[1] - paddleWidth / 2.0;
      topCorners[5][1] = topPaddleCenter[1] + paddleWidth / 2.0;
      topCorners[6][1] = topPaddleCenter[1] + paddleWidth / 2.0;
      topCorners[7][1] = topPaddleCenter[1] - paddleWidth / 2.0;

      topCorners[0][2] = topPaddleCenter[2] - paddleHeight / 2.0;
      topCorners[1][2] = topPaddleCenter[2] - paddleHeight / 2.0;
      topCorners[2][2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
      topCorners[3][2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
      topCorners[4][2] = topPaddleCenter[2] - paddleHeight / 2.0;
      topCorners[5][2] = topPaddleCenter[2] - paddleHeight / 2.0;
      topCorners[6][2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;
      topCorners[7][2] = topPaddleCenter[2] - paddleHeight / 2.0 + paddleLipHeight;

      // rotated and translated
      double topCornersRot[8][3];
      for (int i = 0; i < 8; i++)
      {
        topCorners[i][0] = topCorners[i][0] + topDisp[0];
        double tempTopCorner = topCorners[i][1]; //J.J. - Bug fix- added.
        topCorners[i][1] = topCorners[i][1] * crot - topCorners[i][2] * srot + topDisp[1];
        topCorners[i][2] = tempTopCorner * srot + topCorners[i][2] * crot + topDisp[2];
      }

      double topPadxhat[3] = {1.0, 0.0, 0.0};
      double topPadyhat[3] = {0.0, crot, srot};
      double topPadzhat[3] = {0.0, srot, -crot};

      int ijk[3];
      double pcoords[3];

      // calculate top paddle bounds
      breastCompressed->ComputeStructuredCoordinates(topPos, ijk, pcoords);
      topBound[0] = ijk[0];
      topBound[1] = ijk[0];
      topBound[2] = ijk[1];
      topBound[3] = ijk[1];
      topBound[4] = ijk[2];
      topBound[5] = ijk[2];

      for (int i = 0; i < 8; i++)
      {
        int inVol = breastCompressed->ComputeStructuredCoordinates(topCorners[i], ijk, pcoords);

        // update bounds
        topBound[0] = (ijk[0] < topBound[0]) ? ijk[0] : topBound[0];
        topBound[1] = (ijk[0] > topBound[1]) ? ijk[0] : topBound[1];
        topBound[2] = (ijk[1] < topBound[2]) ? ijk[1] : topBound[2];
        topBound[3] = (ijk[1] > topBound[3]) ? ijk[1] : topBound[3];
        topBound[4] = (ijk[2] < topBound[4]) ? ijk[2] : topBound[4];
        topBound[5] = (ijk[2] > topBound[5]) ? ijk[2] : topBound[5];

        if (!inVol)
        {
          // outside boundary
          topBound[0] = (topCorners[i][0] < breastCompressedBounds[0]) ? breastCompressedExtent[0] : topBound[0];
          topBound[1] = (topCorners[i][0] > breastCompressedBounds[1]) ? breastCompressedExtent[1] : topBound[1];
          topBound[2] = (topCorners[i][1] < breastCompressedBounds[2]) ? breastCompressedExtent[2] : topBound[2];
          topBound[3] = (topCorners[i][1] > breastCompressedBounds[3]) ? breastCompressedExtent[3] : topBound[3];
          topBound[4] = (topCorners[i][2] < breastCompressedBounds[4]) ? breastCompressedExtent[4] : topBound[4];
          topBound[5] = (topCorners[i][2] > breastCompressedBounds[5]) ? breastCompressedExtent[5] : topBound[5];
        }
      }

      // segment top paddle

      // paddle has 3 parts, the main, the transition, and the lip. The
      // transition has a slope on the outside of the paddle.
      // NOTE: because the top paddle's z-direction is negated, it is treated
      // exactly the same as the bottom paddle in this code, and the slope therefore
      // is the same for both.
      double transitionSlope = (paddleLipHeight - paddleHeight) / paddleEdgeRadius;

// If paddle getting cut off in VTI, may need to unparallelize
#pragma omp parallel for
      for (int i = topBound[0]; i <= topBound[1]; i++)
      {
        for (int j = topBound[2]; j <= topBound[3]; j++)
        {
          for (int k = topBound[4]; k <= topBound[5]; k++)
          {
            // position
            int idx[3] = {i, j, k};
            double pos[3];
            breastCompressed->GetPoint(breastCompressed->ComputePointId(idx), pos);

            // position in paddle coordinates
            double pos2[3];
            pos2[0] = pos[0] - topPos[0];
            pos2[1] = (pos[1] - topPos[1]) * topPadyhat[1] + (pos[2] - topPos[2]) * topPadyhat[2];
            pos2[2] = (pos[1] - topPos[1]) * topPadzhat[1] + (pos[2] - topPos[2]) * topPadzhat[2];

            if (fabs(pos2[1]) <= paddleWidth / 2.0)
            {
              // test if in main part of paddle
              if ((fabs(pos2[2]) <= paddleHeight / 2.0) && (pos2[0] < paddleLength / 2.0) && (pos2[0] > -paddleLength / 2.0 + 2.0 * paddleEdgeRadius))
              {
                unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(idx));
                if (*q == tissue.bg)
                {
                  *q = tissue.paddle;
                }
              }
              else if ((pos2[0] < -paddleLength / 2.0 + 2.0 * paddleEdgeRadius) &&
                       (pos2[0] > -paddleLength / 2.0 + paddleEdgeRadius) &&
                       (pos2[2] < paddleHeight / 2.0) &&
                       (pos2[2] > transitionSlope * (pos2[0] - (-paddleLength / 2.0 + paddleEdgeRadius)) + paddleHeight / 2.0 - paddleLipHeight))
              {
                // in transition region
                unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(idx));
                if (*q == tissue.bg)
                {
                  *q = tissue.paddle;
                }
              }
              else if ((pos2[0] < -paddleLength / 2.0 + paddleEdgeRadius) && (pos2[0] > -paddleLength / 2.0) && (pos2[2] > -paddleLipHeight + paddleHeight / 2.0))
              {
                // in lip region

                // test if in curved region
                double dx = -paddleLength / 2.0 + paddleEdgeRadius - pos2[0];
                double pThick = paddleLipHeight - paddleEdgeRadius + sqrt(paddleEdgeRadius * paddleEdgeRadius - dx * dx);
                if (pos2[2] <= pThick + paddleHeight / 2.0 - paddleLipHeight)
                {
                  // in curved region
                  unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(idx));
                  if (*q == tissue.bg)
                  {
                    *q = tissue.paddle;
                  }
                }
              }
            }
          }
        }
      }

      // bounding box for bottom paddle
      int botBound[6];

      double botCorners[8][3];
      // unrotated
      botCorners[0][0] = bottomPaddleCenter[0] - paddleLength / 2.0;
      botCorners[1][0] = bottomPaddleCenter[0] - paddleLength / 2.0;
      botCorners[2][0] = bottomPaddleCenter[0] - paddleLength / 2.0;
      botCorners[3][0] = bottomPaddleCenter[0] - paddleLength / 2.0;
      botCorners[4][0] = bottomPaddleCenter[0] + paddleLength / 2.0;
      botCorners[5][0] = bottomPaddleCenter[0] + paddleLength / 2.0;
      botCorners[6][0] = bottomPaddleCenter[0] + paddleLength / 2.0;
      botCorners[7][0] = bottomPaddleCenter[0] + paddleLength / 2.0;

      botCorners[0][1] = bottomPaddleCenter[1] - paddleWidth / 2.0;
      botCorners[1][1] = bottomPaddleCenter[1] + paddleWidth / 2.0;
      botCorners[2][1] = bottomPaddleCenter[1] + paddleWidth / 2.0;
      botCorners[3][1] = bottomPaddleCenter[1] - paddleWidth / 2.0;
      botCorners[4][1] = bottomPaddleCenter[1] - paddleWidth / 2.0;
      botCorners[5][1] = bottomPaddleCenter[1] + paddleWidth / 2.0;
      botCorners[6][1] = bottomPaddleCenter[1] + paddleWidth / 2.0;
      botCorners[7][1] = bottomPaddleCenter[1] - paddleWidth / 2.0;

      botCorners[0][2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
      botCorners[1][2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
      botCorners[2][2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
      botCorners[3][2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
      botCorners[4][2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
      botCorners[5][2] = bottomPaddleCenter[2] + paddleHeight / 2.0;
      botCorners[6][2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;
      botCorners[7][2] = bottomPaddleCenter[2] + paddleHeight / 2.0 - paddleLipHeight;

      // rotated and translated
      double botCornersRot[8][3];
      for (int i = 0; i < 8; i++)
      {
        botCorners[i][0] = botCorners[i][0] + bottomDisp[0];
        double tempBotCorner = botCorners[i][1];
        botCorners[i][1] = botCorners[i][1] * crot - botCorners[i][2] * srot + bottomDisp[1];
        botCorners[i][2] = tempBotCorner * srot + botCorners[i][2] * crot + bottomDisp[2];
      }

      double botPadxhat[3] = {1.0, 0.0, 0.0};
      double botPadyhat[3] = {0.0, crot, srot};
      double botPadzhat[3] = {0.0, -srot, crot};

      // calculate top paddle bounds
      breastCompressed->ComputeStructuredCoordinates(bottomPos, ijk, pcoords);
      botBound[0] = ijk[0];
      botBound[1] = ijk[0];
      botBound[2] = ijk[1];
      botBound[3] = ijk[1];
      botBound[4] = ijk[2];
      botBound[5] = ijk[2];

      for (int i = 0; i < 8; i++)
      {
        int inVol = breastCompressed->ComputeStructuredCoordinates(botCorners[i], ijk, pcoords);

        // update bounds
        botBound[0] = (ijk[0] < botBound[0]) ? ijk[0] : botBound[0];
        botBound[1] = (ijk[0] > botBound[1]) ? ijk[0] : botBound[1];
        botBound[2] = (ijk[1] < botBound[2]) ? ijk[1] : botBound[2];
        botBound[3] = (ijk[1] > botBound[3]) ? ijk[1] : botBound[3];
        botBound[4] = (ijk[2] < botBound[4]) ? ijk[2] : botBound[4];
        botBound[5] = (ijk[2] > botBound[5]) ? ijk[2] : botBound[5];

        if (!inVol)
        {
          // outside boundary
          botBound[0] = (botCorners[i][0] < breastCompressedBounds[0]) ? breastCompressedExtent[0] : botBound[0];
          botBound[1] = (botCorners[i][0] > breastCompressedBounds[1]) ? breastCompressedExtent[1] : botBound[1];
          botBound[2] = (botCorners[i][1] < breastCompressedBounds[2]) ? breastCompressedExtent[2] : botBound[2];
          botBound[3] = (botCorners[i][1] > breastCompressedBounds[3]) ? breastCompressedExtent[3] : botBound[3];
          botBound[4] = (botCorners[i][2] < breastCompressedBounds[4]) ? breastCompressedExtent[4] : botBound[4];
          botBound[5] = (botCorners[i][2] > breastCompressedBounds[5]) ? breastCompressedExtent[5] : botBound[5];
        }
      }

// segment bottom paddle
// if paddle getting cut off in VTI, may need to unparallelize
#pragma omp parallel for
      for (int i = botBound[0]; i <= botBound[1]; i++)
      {
        for (int j = botBound[2]; j <= botBound[3]; j++)
        {
          for (int k = botBound[4]; k <= botBound[5]; k++)
          {
            // position
            int idx[3] = {i, j, k};
            double pos[3];
            breastCompressed->GetPoint(breastCompressed->ComputePointId(idx), pos);

            // position in paddle coordinates
            double pos2[3];
            pos2[0] = pos[0] - bottomPos[0];
            pos2[1] = (pos[1] - bottomPos[1]) * botPadyhat[1] + (pos[2] - bottomPos[2]) * botPadyhat[2];
            pos2[2] = (pos[1] - bottomPos[1]) * botPadzhat[1] + (pos[2] - bottomPos[2]) * botPadzhat[2];

            if (fabs(pos2[1]) <= paddleWidth / 2.0)
            {
              // test if in main part of paddle
              if ((fabs(pos2[2]) <= paddleHeight / 2.0) && (pos2[0] < paddleLength / 2.0) && (pos2[0] > -paddleLength / 2.0 + 2.0 * paddleEdgeRadius))
              {
                unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(idx));
                if (*q == tissue.bg)
                {
                  *q = tissue.paddle;
                }
              }
              else if ((pos2[0] < -paddleLength / 2.0 + 2.0 * paddleEdgeRadius) &&
                       (pos2[0] > -paddleLength / 2.0 + paddleEdgeRadius) &&
                       (pos2[2] < paddleHeight / 2.0) &&
                       (pos2[2] > transitionSlope * (pos2[0] - (-paddleLength / 2.0 + paddleEdgeRadius)) + paddleHeight / 2.0 - paddleLipHeight))
              {
                // in transition region
                unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(idx));
                if (*q == tissue.bg)
                {
                  *q = tissue.paddle;
                }
              }
              else if ((pos2[0] < -paddleLength / 2.0 + paddleEdgeRadius) && (pos2[0] > -paddleLength / 2.0) && (pos2[2] > -paddleLipHeight + paddleHeight / 2.0))
              {
                // in lip region

                // test if in curved region
                double dx = -paddleLength / 2.0 + paddleEdgeRadius - pos2[0];
                double pThick = paddleLipHeight - paddleEdgeRadius + sqrt(paddleEdgeRadius * paddleEdgeRadius - dx * dx);
                if (pos2[2] <= pThick + paddleHeight / 2.0 - paddleLipHeight)
                {
                  // in curved region
                  unsigned char *q = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(idx));
                  if (*q == tissue.bg)
                  {
                    *q = tissue.paddle;
                  }
                }
              }
            }
          }
        }
      }
    }
    // finished labeling, save compressed breast

    char outputHdrFilename[256];
    char outputLesionFilename[256];
    char outputgzFilename[256];

    if (done)
    {
      // crop the phantom
      int cExt[6];
      int trueExtent[6];
      breastCompressed->GetExtent(cExt);

      for (int i = cExt[0]; i < cExt[1]; i++)
      {
        for (int j = cExt[2]; j < cExt[3]; j++)
        {
          for (int k = cExt[4]; k < cExt[5]; k++)
          {
            unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(i, j, k));
            if (p[0] != tissue.bg && p[0] != tissue.paddle)
            {
              trueExtent[0] = i;
              goto founde0;
            }
          }
        }
      }

    founde0:
      for (int i = cExt[1]; i >= cExt[0]; i--)
      {
        for (int j = cExt[2]; j < cExt[3]; j++)
        {
          for (int k = cExt[4]; k < cExt[5]; k++)
          {
            unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(i, j, k));
            if (p[0] != tissue.bg && p[0] != tissue.paddle)
            {
              trueExtent[1] = i;
              goto founde1;
            }
          }
        }
      }

    founde1:
      for (int j = cExt[2]; j < cExt[3]; j++)
      {
        for (int i = cExt[0]; i < cExt[1]; i++)
        {
          for (int k = cExt[4]; k < cExt[5]; k++)
          {
            unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(i, j, k));
            if (p[0] != tissue.bg && p[0] != tissue.paddle)
            {
              trueExtent[2] = j;
              goto founde2;
            }
          }
        }
      }

    founde2:
      for (int j = cExt[3]; j >= cExt[2]; j--)
      {
        for (int i = cExt[0]; i < cExt[1]; i++)
        {
          for (int k = cExt[4]; k < cExt[5]; k++)
          {
            unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(i, j, k));
            if (p[0] != tissue.bg && p[0] != tissue.paddle)
            {
              trueExtent[3] = j;
              goto founde3;
            }
          }
        }
      }

    founde3:
      for (int k = cExt[4]; k < cExt[5]; k++)
      {
        for (int j = cExt[2]; j < cExt[3]; j++)
        {
          for (int i = cExt[0]; i < cExt[1]; i++)
          {
            unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(i, j, k));
            if (p[0] != tissue.bg && p[0] != tissue.paddle)
            {
              trueExtent[4] = k;
              goto founde4;
            }
          }
        }
      }

    founde4:
      for (int k = cExt[5]; k >= cExt[4]; k--)
      {
        for (int j = cExt[2]; j < cExt[3]; j++)
        {
          for (int i = cExt[0]; i < cExt[1]; i++)
          {
            unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer(i, j, k));
            if (p[0] != tissue.bg && p[0] != tissue.paddle)
            {
              trueExtent[5] = k;
              goto founde5;
            }
          }
        }
      }

    founde5:
      // number of surrounding air voxels to leave on each side
      int buffer = 100;

      int newExtent[6] = {trueExtent[0] - buffer, trueExtent[1] + buffer, trueExtent[2] - buffer,
                          trueExtent[3] + buffer, trueExtent[4] - buffer, trueExtent[5] + buffer};

      for (int i = 0; i < 3; i++)
      {
        if (newExtent[2 * i] < cExt[2 * i])
        {
          newExtent[2 * i] = cExt[2 * i];
        }
        if (newExtent[2 * i + 1] > cExt[2 * i + 1])
        {
          newExtent[2 * i + 1] = cExt[2 * i + 1];
        }
      }

      // crop phantom
      breastCompressed->Crop(newExtent);

      sprintf(outputHdrFilename, "%s/pc_%d.mhd", workDir.c_str(), seed);
      sprintf(outputgzFilename, "%s/pc_%d.raw.gz", workDir.c_str(), seed);
      sprintf(outputLesionFilename, "%s/pc_%d.loc", workDir.c_str(), seed);
    }
    else
    {
      sprintf(outputHdrFilename, "%s/p_%d_%d.mhd", workDir.c_str(), seed, meshCount);
      sprintf(outputgzFilename, "%s/p_%d_%d.raw.gz", workDir.c_str(), seed, meshCount);
      sprintf(outputLesionFilename, "%s/p_%d_%d.loc", workDir.c_str(), seed, meshCount);
    }

    // save locations
    FILE *outLoc;
    outLoc = fopen(outputLesionFilename, "w");
    for (int j = 0; j < newLesionLoc->GetNumberOfPoints(); j++)
    {
      double pt[3];
      newLesionLoc->GetPoint(j, pt);
      fprintf(outLoc, "%E,%E,%E\n", pt[0], pt[1], pt[2]);
    }

    fclose(outLoc);

    // save gz data
    int finalDim[3];
    breastCompressed->GetDimensions(finalDim);

    gzFile gzoutf = gzopen(outputgzFilename, "wb");
    if (gzoutf == NULL)
    {
      cerr << "Unable to open gzip file for writing\n";
    }
    else
    {
      gzbuffer(gzoutf, finalDim[1] * finalDim[2]);
      unsigned char *p = static_cast<unsigned char *>(breastCompressed->GetScalarPointer());
      for (int i = 0; i < finalDim[0]; i++)
      {
        gzwrite(gzoutf, static_cast<const void *>(p), finalDim[1] * finalDim[2]);
        p += finalDim[1] * finalDim[2];
      }

      gzclose(gzoutf);
    }

    // save mhd header
    FILE *hdrFile = fopen(outputHdrFilename, "w");

    if (hdrFile == NULL)
    {
      cerr << "Unable to open mhd file for writing\n";
    }
    else
    {
      fprintf(hdrFile, "ObjectType = Image\n");
      fprintf(hdrFile, "NDims = 3\n");
      fprintf(hdrFile, "BinaryData = True\n");
      fprintf(hdrFile, "BinaryDataByteOrderMSB = False\n");
      fprintf(hdrFile, "CompressedData = False\n");
      fprintf(hdrFile, "TransformMatrix = 1 0 0 0 1 0 0 0 1\n");
      double finalOrigin[3];
      double finalSpacing[3];
      breastCompressed->GetSpacing(finalSpacing);
      // final origin is location of corner voxel because mhd doesn't allow negative voxel indicies
      breastCompressed->GetOrigin(finalOrigin);
      finalOrigin[0] += breastCompressedExtent[0] * finalSpacing[0];
      finalOrigin[1] += breastCompressedExtent[2] * finalSpacing[1];
      finalOrigin[2] += breastCompressedExtent[4] * finalSpacing[2];
      fprintf(hdrFile, "Offset = %6.4f %6.4f %6.4f\n", finalOrigin[0], finalOrigin[1], finalOrigin[2]);
      fprintf(hdrFile, "CenterOfRotation = 0 0 0\n");
      fprintf(hdrFile, "ElementSpacing = %6.4f %6.4f %6.4f\n", finalSpacing[0], finalSpacing[1], finalSpacing[2]);
      fprintf(hdrFile, "DimSize = %d %d %d\n", finalDim[0], finalDim[1], finalDim[2]);
      fprintf(hdrFile, "AnatomicalOrientation = ???\n");
      fprintf(hdrFile, "ElementType = MET_UCHAR\n");
      fprintf(hdrFile, "ObjectType = Image\n");
      if (done)
      {
        fprintf(hdrFile, "ElementDataFile = pc_%d.raw\n", seed);
      }
      else
      {
        fprintf(hdrFile, "ElementDataFile = p_%d_%d.raw\n", seed, meshCount);
      }
      fclose(hdrFile);
    }

    if (autoRemesh)
    {
      std::cout << "################################################################################\n";
      std::cout << "FINISHED MESHING ITERATION " << meshCount + 1 << " of maximum " << MAX_MESHES << "\n";
      std::cout << "################################################################################\n";
    }
    else
    {
      std::cout << "################################################################################\n";
      std::cout << "FINISHED MESHING ITERATION " << meshCount + 1 << " OUT OF " << numMeshes << " \n";
      std::cout << "################################################################################\n";
    }
    meshCount++;

    if (done)
    {
      std::cout << "\nBreast Compression Successful\n";
    }
  }
  return (0);
}
