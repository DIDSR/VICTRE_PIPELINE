/*! \file breastMass.hxx
 *  \brief breastMass header file
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

#ifndef BREASTMASS_HXX_
#define BREASTMASS_HXX_

// boost
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>

#include <omp.h>

// vtk
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPNGWriter.h>
#include <vtkImageData.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkMinimalStandardRandomSequence.h>
#include <vtkBoxMuellerRandomSequence.h>
#include <vtkImageInterpolator.h>

// create spiculation segments
void createBranch(double, double, double, double, double, double, double, double, 
	double, vtkImageData*, vtkMinimalStandardRandomSequence*, vtkBoxMuellerRandomSequence*, boost::program_options::variables_map);


#endif /* BREASTMASS_HXX_ */
