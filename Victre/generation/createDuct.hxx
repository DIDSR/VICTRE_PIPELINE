/*
 * createDuct.hxx
 *
 *  Created on: Dec 22, 2014
 *      Author: cgg
 */

#ifndef CREATEDUCT_HXX_
#define CREATEDUCT_HXX_

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

#ifndef __TISSUESTRUCT__
	#define __TISSUESTRUCT__
	#include "tissueStruct.hxx"
#endif

#ifndef __DUCT__
	#define __DUCT__
	#include "duct.hxx"
#endif

void generate_duct(vtkImageData* breast, boost::program_options::variables_map vm, vtkPoints* TDLUloc, vtkDoubleArray* TDLUattr, 
	unsigned char compartmentId, int* boundBox, tissueStruct* tissue, double* sposPtr, double* sdirPtr, int seed);

#endif /* CREATEDUCT_HXX_ */
