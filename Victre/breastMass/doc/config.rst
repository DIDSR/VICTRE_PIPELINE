Configuration File Parameters
=============================

The configuration file contains all the parameters necessary to generate a breast mass.  They are divided into several categories, primarily by tissue type.
In the following tables each parameter is given a brief description.

base parameters
----------------
============ ========== ======================================================
Name         Type       Notes
============ ========== ======================================================
base.imgRes  float (mm) voxel size
base.seed    integer    random number seed (chosen randomly if not specified)
============ ========== ======================================================

mass parameters
----------------

=================== ============= ====================================================
Name                Type          Notes
=================== ============= ====================================================
mass.lMax           integer       maximum spherical harmonic order
mass.alpha          double (mm)   mean mass radius (mm)")
mass.meanSigma2     double (mm^2) mean mass surface irregularity variance
mass.stdSigma2      double (mm^2) mass surface irregularity standard deviation
mass.powerLaw       double        covariance power law index
mass.meanLF         double        mean number of low freq. modifications
mass.stdLF          double        std. deviation of number of low freq. modifications
mass.meanShape      double        mean LF shape distribution 0=spike, 1=bump
mass.stdShape       double        std. deviation of LF shape 0=spike, 1=bump
mass.meanLFRad      double        mean relative LF radius
mass.stdLFRad       double        std. deviation of LF radius
mass.meanLFLen      double        mean relative LF length
mass.stdLFLen       double        std. deviation of LF length
mass.meanFuzzAlpha  double        mean fuzzy alpha
mass.stdFuzzAlpha   double        std. deviation of fuzzy alpha
=================== ============= ====================================================

spicule parameters
------------------

====================== ======= ===============================================
Name                   Type    Notes
====================== ======= ===============================================
spicule.meanInitial    double  mean number of initial segments
spicule.stdInitial     double  std. deviation number of initial segments
spicule.meanNeigh      double  mean max number of neighbor segments
spicule.stdNeigh       double  std. deviation max number of neighbor segments
spicule.meanInitRad    double  mean initial relative radius
spicule.stdInitRad     double  std. deviation initial relative radius
spicule.meanRadDec     double  mean radius decrease
spicule.stdRadDec      double  std. deviation radius decrease
spicule.meanInitLen    double  mean initial relative radius
spicule.stdInitLen     double  std. deviation initial relative radius
spicule.meanLenDec     double  mean length decrease
spicule.stdLenDec      double  std. deviation length decrease
spicule.meanContProb   double  mean continue prob.
spicule.stdContProb    double  std. deviation continue prob.
spicule.meanBranchAng  double  mean asymmentric bifurcation prob.
spicule.stdBranchAng   double  std. deviation asymm bif. prob.
====================== ======= ===============================================
