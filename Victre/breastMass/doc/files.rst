Output Files
============

breastMass outputs a 3-dimensional voxel volume consisting of binary values (1 for inside mass and 0 for outside mass).  Three output files are generated:

mass\_\ *nnnnnnnn*.vti
   The voxel volume stored in VTK image format (XML with zlib-compressed data element).  This is a file format native to the Visualization Toolkit.  See https://www.vtk.org for
   more information.  This file can be opened in Paraview, an open-source, multi-platform data analysis and visualization application available for download at
   https://www.paraview.org/  The integer *nnnnnnnn* appearing in the filename is the random number seed that was used to generate the phantom.

mass\_\ *nnnnnnnn*\_\ *mmm*.raw
   The voxel volume stored in raw binary format.  The integer *nnnnnnnn* appearing in the filename is the random number seed that was used to generate the phantom.
   The integer *mmm* is the number of voxels in each dimension (all voxel volumes are cubes).
   
mass\_\ *nnnnnnnn*.cfg
   A copy of the configuration file used to generate the mass.

