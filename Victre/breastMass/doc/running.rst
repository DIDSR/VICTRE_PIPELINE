Running the Software
====================

Command line options
--------------------

breastMass is run from the command line and requires that a single command line options be specified::

  > breastMass -c <configFile>

required command line argument

-c configFile       a configuration file specifying all mass parameters

optional command line argument

-h                  print a help message

Operating Mode
--------------

breastMass is run from the command line.  A configuration file must be provided as in the following example::

  > breastCrop -c mass.cfg

See subsequent sections of this documentation for the output file formats and configuration file details.

Sample Configuration Files
--------------------------

The /cfg directory contains several sample configuration files

**testMass.cfg**
  A random 20mm diameter mass without spicules generated at 0.1 mm resolution

**testSpiculated.cfg**
  A random 20mm diameter mass with spicules generated at 0.1 mm resolution

**VICTREMass.cfg**
  The 5mm diameter spiculed mass used in the VICTRE study with specific random number seed (0.05 mm voxel resolution)



