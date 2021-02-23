/*! \file breastCompress.hxx
 *  \brief breastCompress main header file
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

#ifndef BREASTCOMPRESS_HXX_
#define BREASTCOMPRESS_HXX_

#include <iostream>
#include <fstream>
#include <cstdint>
#include <math.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>

#include <zlib.h>

#include <boost/program_options.hpp>

#ifndef __OMP__
#define __OMP__
#include <omp.h>
#endif

// vtk stuff
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkMath.h>
#include <vtkPolygon.h>
#include <vtkPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkGeometryFilter.h>
#include <vtkContourFilter.h>
#include <vtkReverseSense.h>
#include <vtkQuadricDecimation.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkCellArray.h>
#include <vtkPointData.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkSurfaceReconstructionFilter.h>
#include <vtkPLYWriter.h>
#include <vtkImageData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkXMLImageDataReader.h>
#include <vtkMetaImageReader.h>
#include <vtkMetaImageWriter.h>
#include <vtkVector.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkDiscreteMarchingCubes.h>
#include <vtkTetra.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkPoints.h>

#include "tetgen.h"

#include "tissueStruct.hxx"

#endif /* BREASTCOMPRESS_HXX_ */
