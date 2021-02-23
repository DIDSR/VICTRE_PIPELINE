Installation
============

The breastCompress source code is available `here <https://github.com/DIDSR/breastCompress>`_.  It is written in C++ and therefore must be compiled for your system before use.

Operating System
----------------

This software has been developed and tested on Debian and Red Hat based Linux distributions.  Compilation is via cmake and should be possible on any UNIX-like operationg system (including MacOS) with minimal
modification to the build configuration in the provided CMakeLists.txt file, however the following distributions are recommended: 

- Ubuntu 16.04
- Centos 7.x

The following build instructions are based on Ubuntu 16.04.

Prerequisites
-------------

* `FEBiO <https://febio.org/febio/>`_ (required)

  To help perform the breast compression simulation this software relies upon elements of the `FEBiO Software Suite <https://febio.org>`_, a nonlinear finite element analysis software package.
  FEBiO is not distributed as part of breastCompress and must be downloaded from the `FEBiO website <https://febio.org/febio/>`_.  Registration is required but
  FEBiO is open source software and free for non-commerical use.  

* `PostView <https://febio.org/postview/>`_ (optional)

  PostView is another component of the FEBiO Sofware Suite which can be used to visualize the compression process.  This is not required to use breastCompress, but may be useful for debugging compression
  problems or producing graphics, strain maps, displacement maps, etc.

Compilation
-----------

First, install the packages necessary to compile the software (reminder: instructions are specific to Ubuntu,
consult your operating system documentation for analagous package names and installation procedures)  From a terminal command prompt::

    > sudo apt-get install cmake vtk6 libvtk6-dev \
			   libblas-dev liblapack-dev \
			   libboost-dev libboost-program-options-dev \
			   libproj-dev libtet1.5-dev zlib1g-dev

Next, clone the git repository::

    > git clone https://github.com/DIDSR/breastCompress

From a suitable build directory generate a Makefile using cmake and run make::

    > cmake <path to repository>
    > make

where <path to repository> is the directory containing the cloned repository.
.. note:: When running cmake you may recieve an error message related to references to files "/usr/lib/x86_64-linux-gnu/libvtkRenderingPythonTkWidgets.so" and "/usr/bin/vtk".
	  This is a known Ubuntu VTK bug but since we are not linking to these files, the error can safely be ignored.

FeBiO Setup
-----------

To function properly breastCompress needs to know where the FeBiO executable is installed on your system.  This information can be passed to FeBiO as a command line option
(see next section) or be specified in the environment variable $FEBIO_PATH.  If using the bash shell this variable can be set using::

  > export FEBIO_PATH=<path>

for example::

  > export FEBIO_PATH=/usr/local/febio-2.3.1

breastCompress will expect the FeBiO executable to be located at $FEBIO_PATH/bin/febio2.lnx64.


Finish
------

If the previous steps were successful, the executable breastPhantom should be present in the build directory.  This file can be copied to a convenient location in your path,
or the build directory can be added to $PATH.  No additional installation steps are required.  Consult the next section which describes how to run the executable.
