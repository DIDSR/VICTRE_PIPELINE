# mass configuration file

[base]
# voxel size (mm)
imgRes=0.05
# complexity scaling for VICTRE 0=simple, 1=full complexity
complexity=1.0
# random number seed (unsigned int) from /dev/urandom if unspecified
seed=$seed

[mass]
# maximum spherical harmonic order
lMax=4
# mean mass radius (mm)
alpha=$alpha
# mean mass surface irregularity variance (mm^2)
meanSigma2=0.31
# mass surface irregularity standard deviation (mm^2)
stdSigma2=0.04
# covariance power law index
powerLaw=4.0
# mean number of low freq. modifications
meanLF=611.2
# std. deviation of number of low freq. modifications
stdLF=70.6
# mean LF shape distribution 0=spike, 1=bump
meanShape=0.36
# std. deviation of LF shape 0=spike, 1=bump
stdShape=0.48
# mean relative LF radius
meanLFRad=0.229
# std. deviation of LF radius
stdLFRad=0.073
# mean relative LF length
meanLFLen=0.113
# std. deviation of LF length
stdLFLen=0.021
# mean fuzzy alpha
meanFuzzAlpha=0.015
# std. deviation of fuzzy alpha
stdFuzzAlpha=0.0

[spicule]
# mean number of initial segments
meanInitial=1358
# std. deviation number of initial segments
stdInitial=365
# mean max number of neighbor segments
meanNeigh=8.98
# std. deviation max number of neighbor segments
stdNeigh=1.89
# mean initial relative radius
meanInitRad=0.024
# std. deviation initial relative radius
stdInitRad=0.0053
# mean radius decrease
meanRadDec=0.89
# std. deviation radius decrease
stdRadDec=0.31
# mean initial relative radius
meanInitLen=0.173
# std. deviation initial relative radius
stdInitLen=0.018
# mean length decrease
meanLenDec=0.91
# std. deviation length decrease
stdLenDec=0.30
# mean continue prob.
meanContProb=0.717
# std. deviation continue prob.
stdContProb=0.057
# mean branch angle.
meanBranchAng=6.55
# std. deviation branch angle.
stdBranchAng=0.62
