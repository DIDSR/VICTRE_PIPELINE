Getting started
===============

Requirements
------------

-  For phantom generation
-  cmake 3+
-  gcc/g++
-  VTK library
-  Lapack library
-  boost library
-  For phantom compression
-  FEBio
-  gcc/g++
-  cmake 3+
-  VTK library
-  Lapack library
-  boost library
-  For mass generation
-  gcc/g++
-  cmake 3+
-  VTK library
-  Lapack library
-  boost library
-  For MCGPU projection
-  CUDA 10.4+
-  NVIDIA GPU (8GB+ recommended)
-  gzip
-  (Optional) openMPI
-  For reconstruction
-  gcc/g++
-  For the unified pipeline script
-  Python 3.6+

Installation
------------

Before proceeding, make sure you have all the requirements listed above.
You can also install the required libraries using this:

``sudo apt-get install cmake vtk7 libvtk7-dev        libblas-dev liblapack-dev libopenmpi-dev       libboost-dev libboost-program-options-dev        libproj-dev zlib1g-dev gzip``

Clone the repository

``git clone https://github.com/DIDSR/VICTRE_PIPELINE``

Execute the installation script:

``source install.sh``

Follow the instructions and compile the 5 parts. You might need to edit
the ``install.sh`` file to change the route of some libraries at the
beginning.

Install the required python libraries:

``pip install numpy scipy termcolor progressbar2 h5py pydicom``

    You might need to use ``pip3`` instead of ``pip``.

Add the ``FEBio`` executable to the path (replace ``{{{routetofebio}}}``
with the path to your FEBio installation):

``export PATH="$PATH:/{{{routetofebio}}}/FEBio-2.9.1/bin"``

    You can also add that line to your ``.bashrc`` file to make it
    permanent. Note the version number for FEBio on the path, it might 
    be different.

Usage
-----

Copy the examples files to the parent folder and run them in a
GPU-enabled machine with CUDA:

``python example1.py``

    You might need to use ``python3`` instead of ``python``.

Each example file starts the pipeline from a different step, use
``example1.py`` to start from the phantom generation (it will need hours
to complete), use ``example5.py`` to start from the projection step (it
will finish in about ~10 minutes depending on your computer). When
finished, you will find the results in the ``results`` folder under the
``1`` subfolder.

-  ``p_1.raw.gz``: original phantom
-  ``pc_1.raw.gz``: compressed phantom
-  ``pc_1_crop.raw.gz``: cropped compressed phantom
-  ``pcl_1.raw.gz``: compressed original phantom with the inserted
   lesions
-  ``pcl_1.loc``: file containing the coordinates of the inserted
   lesions in the phantom
-  Last number is the lesion type: ``1`` for calcification clusters,
   ``2`` for masses
-  ``projection_DM1.raw``: contains the DM projection with a size of
   3000 x 1000 pixels in raw format
-  ``reconstruction1.raw``: contains the DBT reconstruction with a size
   of 3000 x 1000 x 25 voxels in raw format
-  ``ROIs.h5``: contains the lesion-present and lesion-absent regions of
   interest.
-  ``ROIs``: subfolder will also contain the ROIs in raw format (size is
   specified in the code, ``109 x 109 x 9`` in the examples)
-  ``ROI_DM_XX_typeT``: DM cropped image for lesion number ``XX`` of
   lesion type ``T`` (absent regions will have ``T < 0``)
-  ``ROI_DBT_XX_typeT``: DBT cropped volume for lesion number ``XX`` of
   lesion type ``T`` (absent regions will have ``T < 0``)
-  ``T = 1`` for calcification clusters, ``T = 2`` for masses

