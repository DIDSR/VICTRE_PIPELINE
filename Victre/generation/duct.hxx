/*
 * duct.hxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#ifndef DUCT_HXX_
#define DUCT_HXX_

#ifndef __CMATH__
#define __CMATH__
#include <cmath>
#endif

#ifndef __OMP__
#define __OMP__
#include <omp.h>
#endif

#ifndef __ALGORITHM__
#define __ALGORITHM__
#include <algorithm>
#endif

#ifndef __VTKIMAGEDATA__
#define __VTKIMAGEDATA__
#include <vtkImageData.h>
#endif

#ifndef __VTKPOINTS__
#define __VTKPOINTS__
#include <vtkPoints.h>
#endif

#ifndef __VTKDOUBLEARRAY__
#define __VTKDOUBLEARRAY__
#include <vtkDoubleArray.h>
#endif

#ifndef __VTKMATH__
#define __VTKMATH__
#include <vtkMath.h>
#endif

#ifndef __VTKVECTOR__
#define __VTKVECTOR__
#include <vtkVector.h>
#endif

#ifndef __BOOST__
#define __BOOST__
#include <boost/random.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/program_options.hpp>
#endif

#ifndef __TISSUESTRUCT__
#define __TISSUESTRUCT__
#include "tissueStruct.hxx"
#endif

// forward declaration
class ductSeg;
class ductBr;

/**********************************************
*
* structure for duct tree initialization
*
**********************************************/
struct ductTreeInit{
  // random number generator seed
  int seed;
  // pointer to bound box
  int *boundBox;
  // compartment id
  unsigned char compartmentId;
  // segmentation tissue values
  tissueStruct* tissue;
  // FOV
  double startPos[3];
  double endPos[3];
  // preferential direction of growth (same as startDir)
  double prefDir[3];
  // size of arrays
  unsigned int nVox[3];
  unsigned int nFill[3];
  // pointer to breast
  vtkImageData* breast;
  // pointer to TDLU locations
  vtkPoints* TDLUloc;
  // pointer to TDLU attributes
  vtkDoubleArray* TDLUattr;
};


/**********************************************
*
* Class for a duct tree
*
**********************************************/

class ductTree{

  friend class ductBr;
  friend class ductSeg;

  typedef boost::mt19937 rgenType;
  // random number generator - constructor should set seed!!
  rgenType randGen;
  // pointer to configuration
  boost::program_options::variables_map opt;
public:
  // pointer to breast bound box
  int *boundBox;
  // compartment id
  unsigned char compartmentId;
  // tissue values
  tissueStruct* tissue;
  // maximum number of branches
  unsigned int maxBranch;
  // base length of initial branch
  double baseLength;
  // fill map giving distance to tree in roi
  // initial value is distance to base of tree
  vtkImageData* fill;
  // duct tree count
  static unsigned int num;
  // uniform [0,1) distribution
  boost::uniform_01<rgenType> u01;
  // beta distribution for segment length
  //boost::math::beta_distribution<> lengthDist;
  // beta distribution for radius of curvature
  boost::math::beta_distribution<> radiusDist;
  // duct tree id number
  unsigned int id;
  // keep track of number of branches in tree
  unsigned int numBranch;
  // pointer to main branch
  ductBr* head;
  // pointer to breast
  vtkImageData* breast;
  // pointer to TDLU locations
  vtkPoints* TDLUloc;
  // pointer to TDLU attributes
  vtkDoubleArray* TDLUattr;
  // preferential growth direction
  double prefDir[3];
  // save to file function
  // constructor
  ductTree(boost::program_options::variables_map, ductTreeInit*);
  // destructor
  ~ductTree();
};



/**********************************************
*
* Class for a duct branch
*
**********************************************/

class ductBr {
  // this is one branch of a tree
  // sType refers to segment type (equivalent to tree type)

  friend class ductSeg;
  friend class ductTree;

  // start and end position of branch
  double startPos[3];
  double endPos[3];
  // start and end radius (mm)
  double startRad, endRad;
  // start and end direction (unit vector)
  double startDir[3];
  double endDir[3];
  // rotation angle from parent
  double azimuth;
  // length of branch and current length (mm)
  double length, curLength;
  // pointer to first segment of branch
  ductSeg* firstSeg;
  // pointer to last segment of branch
  ductSeg* lastSeg;
  // pointer to parent branch
  ductBr* parent;
  // pointer to first child branch
  ductBr* firstChild;
  // pointer to second child branch
  ductBr* secondChild;
  // pointer to sibling branch
  ductBr* sibBranch;
  // branch id number
  unsigned int id;
  // number of child branches (0 or 2)
  unsigned int nChild;
  // pointer to tree instance
  ductTree* myTree;
  // level in network, 0 == main branch
  unsigned int level;
  // generation of branch, 0 == root
  unsigned int gen;
  // function to set length of branch
  double setLength(void);
  // function to set number of children
  unsigned int setChild(void);
  // function to pick starting radii of child branches
  void setRadiiThetas(double*,double*);
  // function to pick starting direction based on parent direction
  void setDir(double*,double);
public:
  // constructor for first branch (the root)
  ductBr(double*, double*, double, ductTree*);
  // constructor for first child branch of a parent branch
  ductBr(ductBr*, unsigned int, unsigned int, double, double);
  // constructor for second child branch
  ductBr(ductBr*, ductBr*, unsigned int, unsigned int, double, double);
  // destructor that deletes all child branches as well
  ~ductBr();
};



/**********************************************
*
* Class for a segment (of a duct branch)
*
**********************************************/

class ductSeg {
  // this is one segment of a branch
  
  friend class ductBr;
  friend class ductTree;

  // start and end position of segment
  double startPos[3];
  double endPos[3];
  // start and end radius rate of change
  double startDeriv, endDeriv;
  // center of curvature
  double centerCurv[3];
  // start and end direction (unit vector)
  double startDir[3];
  double endDir[3];
  // radius of curvature
  double radCurv;
  // pointer to previous segment of branch
  ductSeg* prevSeg;
  // pointer to owning branch
  ductBr* myBranch;
  // cubic spline coefficients
  double shape[4];
public:
  // start and end radius (mm)
  double startRad, endRad;
  // length of segment (mm)
  double length;
  // pointer to next segment of branch
  ductSeg* nextSeg;
  // make a first segment - determines endPos, endRad
  // centerCurv, length, endDir
  void makeSeg(void);
  // make shape
  void setShape(void);
  // get segment radius
  double getRadius(double);
  // update voxel-based map of duct tree - this edits breast data
  void updateMap(void);
  // constructor for first segment
  ductSeg(ductBr*);
  // constructor for subsequent segments
  ductSeg(ductSeg*);
};





#endif /* DUCT_HXX_ */
