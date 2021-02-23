/*
 * vein_options.hxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#ifndef VEIN_HXX_
#define VEIN_HXX_

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

#ifndef __VTKMATH__
#define __VTKMATH__
#include <vtkMath.h>
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
class veinSeg;
class veinBr;

/**********************************************
*
* structure for veinTree initialization
*
**********************************************/
struct veinTreeInit{
  // random number generator seed
  int seed;
  // pointer to bound box
  int *boundBox;
  // tissue values
  tissueStruct* tissue;
  // FOV
  double startPos[3];
  double endPos[3];
  // preferential direction of growth (same as startDir)
  double nipplePos[3];
  // size of arrays
  unsigned int nVox[3];
  unsigned int nFill[3];
  // pointer to breast
  vtkImageData* breast;
};


/**********************************************
*
* Class for an vein tree
*
**********************************************/

class veinTree {

  friend class veinBr;
  friend class veinSeg;

  typedef boost::mt19937 rgenType;
  // random number generator - constructor should set seed!!
  rgenType randGen;
  // pointer to configuration
  boost::program_options::variables_map opt;
public:
  // pointer to breast bound box
  int *boundBox;
  // tissue values
  tissueStruct* tissue;
  // maximum number of branches
  unsigned int maxBranch;
  // length of base branch
  double baseLength;
  // fill map giving distance to tree in roi
  // initial value is distance to base of tree
  vtkImageData* fill;
  // vein tree count
  static unsigned int num;
  // uniform [0,1) distribution
  boost::uniform_01<rgenType> u01;
  // beta distribution for radius of curvature
  boost::math::beta_distribution<> radiusDist;
  // vein tree id number
  unsigned int id;
  // keep track of number of branches in tree
  unsigned int numBranch;
  // pointer to main branch
  veinBr* head;
  // pointer to breast
  vtkImageData* breast;
  // preferential growth direction
  double nipplePos[3];
  // save to file function
  // constructor
  veinTree(boost::program_options::variables_map, veinTreeInit*);
  // destructor
  ~veinTree();
};



/**********************************************
*
* Class for a vein branch
*
**********************************************/

class veinBr {
  // this is one branch of a tree
  // sType refers to segment type (equivalent to tree type)

  friend class veinSeg;
  friend class veinTree;

  // start and end position of branch
  double startPos[3];
  double endPos[3];
  // start and end radius (mm)
  double startRad, endRad;
  // start and end direction (unit vector)
  double startDir[3];
  double endDir[3];
  // length of branch and current length (mm)
  double length, curLength;
  // rotation angle from parent
  double azimuth;
  // pointer to first segment of branch
  veinSeg* firstSeg;
  // pointer to last segment of branch
  veinSeg* lastSeg;
  // pointer to parent branch
  veinBr* parent;
  // pointer to first child branch
  veinBr* firstChild;
  // pointer to second child branch
  veinBr* secondChild;
  // pointer to sibling branch
  veinBr* sibBranch;
  // branch id number
  unsigned int id;
  // number of child branches
  unsigned int nChild;
  // pointer to tree instance
  veinTree* myTree;
  // level in network, 0 == main branch
  unsigned int level;
  // generation of branch, 0 == root
  unsigned int gen;
  // function to set length of branch
  double setLength(void);
  // function to set number of children
  unsigned int setChild(void);
  // function to pick starting radii and angles of child branches
  void setRadiiThetas(double*,double*);
  // function to pick starting direction based on parent direction
  void setDir(double*, double);
public:
  // constructor for first branch (the root)
  veinBr(double*, double*, double, veinTree*);
  // constructor for first child branch of a parent branch
  veinBr(veinBr*, unsigned int, unsigned int, double, double);
  // constructor for other branches
  veinBr(veinBr*, veinBr*, unsigned int, unsigned int, double, double);
  // destructor that deletes all child branches as well
  ~veinBr();
};



/**********************************************
*
* Class for a segment (of a vein branch)
*
**********************************************/

class veinSeg {
  // this is one segment of a branch

  friend class veinBr;
  friend class veinTree;

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
  veinSeg* prevSeg;
  // pointer to owning branch
  veinBr* myBranch;
  // cubic spline coefficients
  double shape[4];
public:
  // start and end radius (mm)
  double startRad, endRad;
  // length of segment (mm)
  double length;
  // pointer to next segment of branch
  veinSeg* nextSeg;
  // make a first segment - determines endPos, endRad
  // centerCurv, length, endDir
  void makeSeg(void);
  // make shape
  void setShape(void);
  // get segment radius
  double getRadius(double);
  // update voxel-based map of vein tree - this edits breast data
  void updateMap(void);
  // constructor for first segment
  veinSeg(veinBr*);
  // constructor for subsequent segments
  veinSeg(veinSeg*);
};

#endif /* VEIN_HXX_ */
